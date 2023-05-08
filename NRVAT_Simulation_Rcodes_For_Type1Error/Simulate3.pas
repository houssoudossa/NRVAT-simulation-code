program simulate3(output,pedfile);
{
 This program version uses an updated random number generator.
 Modified SIMULATE - for SNPs, allows autosomal SNPs to be in LD.

 Simulation program. Output from this program is in SLINK format, and
 is ready for analysis with ISIM, LSIM, or MSIM of the SLINK package.

 The files "problem" through "simout" are defined locally; see below.
 Only "pedfile" is defined globally.

 File     External name   i/o  Defined in routine
 ------------------------------------------------
 pedfile   pedfile.dat    write    main
 problem   problem.dat    read     readfile
 problem   problem.dat    write    writefile
 simdata   simdata.dat    read     readdata
 simped    simped.dat     read     readped
 simout    simout.dat     write    writesim
}

type
 real=single;
 integer=longint;

label 99;

const
 version='8 Sep 2011';
 maxpeds=10;		{Max. no. of pedigrees per replicate}
 maxloci=2000;		{Max. number of loci}
 maxind=50;			{Max. no. individuals in one pedigree}
 maxall=2;			{Maximum number of alleles at any one locus; do not change}
 maxgen = (maxall*(maxall+1)) div 2;
 maxfact=maxall;
 verbose=true;

var
 low : 1..2;
 pedfile:text;
 gasdeviset,sexlink : 0..1;
 sexdiff: 0..2;
 sexratio:real;
 gasdevgset:real;
 doubled:array[1..maxpeds,1..5,1..2] of 0..maxind;
 numgen,numreps,ss,numpeds,numloci:integer;
 liabnum,locustype,numall:array[1..maxloci] of integer;
 theta:array[0..maxloci,1..2] of real;
 genefreq:array[1..maxloci,1..maxall] of real;
 numind: array[1..maxpeds] of integer;
 indnum,pa,ma,fo,nps,nms,pro,disease,
  liabcl:array[1..maxpeds,1..maxind] of integer;
 sex:array[1..maxpeds,1..maxind] of 1..2;
 typed:array[1..maxpeds,1..maxind,1..maxloci] of 0..1;
 locus:array[1..maxpeds,1..maxind,1..maxloci,1..2] of 1..maxall;
 simmed:array[1..maxpeds,0..maxind] of boolean;
 nfact,ntrait:array[1..maxloci] of integer;
 multvar,vari:array[1..maxloci] of real;
 mean,pen:array[1..maxloci,1..maxgen] of real;
 binall:array[1..maxloci,1..maxall,1..maxfact] of 0..1;
 facttot:array[1..maxloci,1..maxfact] of 0..1;
 qphen:array[1..maxpeds,1..maxind,1..maxloci] of real;
 aphen:array[1..maxpeds,1..maxind,1..maxloci] of -1..2;
 bad:boolean;
 seqrun:integer;
 corr2:real; {r^2 for LD between consecutive SNPs}

 {For random number generator}
 u,v,w:QWord;
 seed:QWord; {Called j in Press book}

{$I prog.p}



PROCEDURE readfile;
VAR
 problem:text;
BEGIN
 if verbose then writeln('Opening "problem.dat" file for input of number of replicates');
 assign(problem,'problem.dat');
 reset(problem);
 read(problem,numreps);
 if (not seekeoln(problem)) then read(problem,seqrun);
 close(problem);
END; {readfile}


PROCEDURE writefile;
{Write number of replicates to PROBLEM.DAT for next time}
VAR
 problem:text;
BEGIN
 if verbose then writeln('Opening "problem.dat" file for output');
 assign(problem,'problem.dat');
 rewrite(problem);
 write(problem,'  ',numreps);
 if seqrun>-1 then begin
  seqrun:=seqrun+1;
  write(problem,'  ',seqrun);
 end;
 writeln(problem);
{
 writeln;
 writeln('Output files created (analogous to SLINK files):');
 writeln('  PEDFILE.DAT   simulated pedigree data');
 writeln('  SIMOUT.DAT    summary on simulation');
 writeln('  PROBLEM.DAT   containing number of replicates');
}
 close(problem);
END; {writefile}




FUNCTION nordev(mu,sig:real):real;

 FUNCTION gasdev:real;
 VAR
  fac,r,v1,v2:real;
 BEGIN
  IF gasdeviset=0 THEN BEGIN
   REPEAT
    v1:=2.0*doub-1.0;
    v2:=2.0*doub-1.0;
    r:=sqr(v1)+sqr(v2);
   UNTIL (r<1.0) AND (r>0.0);
   fac:=sqrt(-2.0*ln(r)/r);
   gasdevgset:=v1*fac;
   gasdev:=v2*fac;
   gasdeviset:=1;
  END ELSE BEGIN
   gasdeviset:=0;
   gasdev:=gasdevgset;
  END;
 END; {gasdev}

VAR
 tmp:real;

BEGIN
 tmp:=sig*gasdev+mu;
 IF abs(tmp)<0.00001 THEN tmp:=0.000001;
 nordev:=tmp;
END; {nordev}


PROCEDURE readdata;
VAR     (* File Must Be Named Simdata.dat *)
 simdata:text;
 i,j,ll,k,loopend:integer; {changed - loopend added}
 totfreq:real; {changed}
BEGIN
 if verbose then writeln('Opening "simdata.dat" file for input');
 assign(simdata,'simdata.dat');
 reset(simdata);
 readln(simdata,numloci,i,sexlink);
 readln(simdata); readln(simdata); (* Loci Must Be In Chromosome Order *)
 FOR i:=1 TO numloci DO BEGIN
  readln(simdata,locustype[i],numall[i]);
  FOR j:=1 TO numall[i] DO
   read(simdata,genefreq[i,j]);
  readln(simdata);

  totfreq:=0;  {Test for sum of gene frequencies = 1.  Changed}
  for j:=1 to numall[i] do totfreq:=totfreq+genefreq[i,j];
  if abs(totfreq-1.0)>0.01 then begin
   writeln('ERROR:  Allele frequencies at locus ',i:1,' do not sum to 1.');
   writeln('        Press Ctrl-C to abort program.');
   readln;
  end;

  CASE locustype[i] OF
   (* First Locus Must Be Affection Status, Other Affection Status
      Loci can have only 1 Liability Class *)
   1:BEGIN
      readln(simdata,liabnum[i]);
      IF sexlink=1 THEN loopend:=liabnum[i]*2 else loopend:=liabnum[i]; {changed}
      FOR j:=1 TO loopend DO IF i>1 THEN BEGIN {changed}
       numgen:=((numall[i])*(numall[i]+1)) div 2; {changed}
       IF sexlink=1 THEN numgen:=numgen+numall[i];
       FOR ll:=1 TO numgen DO read(simdata,pen[i,ll]);
       readln(simdata);
      END
       ELSE readln(simdata);
     END;
   2:BEGIN
      readln(simdata,nfact[i]);
      FOR j:=1 TO numall[i] DO
       FOR k:=1 TO nfact[i] DO
	read(simdata,binall[i,j,k]);
      readln(simdata);
     END;
   0:BEGIN
      readln(simdata,ntrait[i]);
      IF ntrait[i]>1 THEN writeln('ONLY ONE QUANTITATIVE TRAIT AT A LOCUS ALLOWED!');
      numgen:=(numall[i]*(numall[i]+1)) div 2;
      FOR j:=1 TO numgen DO
       read(simdata,mean[i,j]);
      readln(simdata);
      readln(simdata,vari[i]);
      readln(simdata,multvar[i]);
     END;
   3:BEGIN
     END;
  END;
 END;
 readln(simdata,sexdiff);
 FOR i:=1 TO (numloci-1) do BEGIN
  read(simdata,theta[i,1]);
  theta[i,2]:=theta[i,1];
 END;
 readln(simdata);
 if sexdiff = 1 then BEGIN
   readln(simdata,sexratio);
   for i:= 1 to (numloci-1) DO
     theta[i,2]:= 0.5*(1-exp(sexratio*ln(1-2*theta[i,1]))); {changed - fewer parenth}
 END ELSE
 if sexdiff = 2 then
  for i:= 1 to (numloci-1) DO
   read(simdata,theta[i,2]);
 readln(simdata);
 close(simdata);
END; {readdata}


PROCEDURE readped(VAR bad:boolean); (* If Disease is First Locus,
	    Given Typings will be used for Affection Status Type *)
LABEL 90;
VAR     (* File Must Be Called Simped.Dat and Must Be in Post-Makeped Format *)
 simped:text; (* File Must Have Header Line With Number of Pedigrees *)
 i,j,k,pedn:integer; (* Followed By Number of People In Each Pedigree *)

BEGIN   (* For Each Locus Besides The Disease, 0 if untyped, 1 if typed *)
 if verbose then writeln('Opening "simped.dat" file for input');
 assign(simped,'simped.dat');
 reset(simped);
 read(simped,numpeds);
 if numpeds>maxpeds then begin
  writeln(' ERROR: Number of pedigrees specified too large (>',maxpeds,')');
  bad:=true;
  goto 90;
 end;
 FOR i:=1 TO numpeds DO read(simped,numind[i]);

 if not seekeoln(simped) then begin
  read(simped,corr2);  {read r^2 if present}
  {writeln('Input value of r^2 =',corr2:7:4,' encountered on line 1 in "simped.dat" file');}
 end;
 readln(simped);

 if numind[1]=1 then begin
  writeln('WARNING: Pedigree 1 contains only 1 individual. Perhaps header line in');
  writeln('         SIMPED.DAT missing. Press <ENTER> to continue or Ctrl-C to stop.');
  readln;
 end;
 low:=1;
 for j:= 1 to maxpeds do
   for i:= 1 to 5 do
      for k:= 1 to 2 do
	 doubled[j,i,k]:=0;
 FOR j:=1 TO numpeds DO
  FOR i:=1 TO numind[j] DO BEGIN
   if numind[j]>maxind then begin
    writeln(' ERROR: number of individuals larger than max of ',maxind:1);
    writeln(j,'-th pedigree');
    bad:=true;
    goto 90;
   end;
   read(simped,pedn,indnum[j,i],pa[j,i],ma[j,i],fo[j,i],nps[j,i],nms[j,i],sex[j,i],pro[j,i]);
   if indnum[j,i] <> i THEN BEGIN
     writeln(' ERROR:  INDIVIDUALS MUST BE NUMBERED SEQUENTIALLY IN EACH PEDIGREE,');
     writeln(' STARTING FROM 1!  Pedigree ',j:1,'  ',pedn:1,'  (top line incorrect?)');
     bad:=true;
     goto 90;
   END;
   if pro[j,i] > 1 then if pa[j,i]=0 then doubled[j,pro[j,i],1]:=i else doubled[j,pro[j,i],2]:=i;
   IF locustype[1]=1 THEN BEGIN read(simped,disease[j,i]);
    IF liabnum[1]>1 THEN read(simped,liabcl[j,i]);
    low:=2;
   END;

   FOR k:=low TO numloci DO read(simped,typed[j,i,k]);
   readln(simped);

  END;

  FOR i:= 1 to numpeds do
   IF (doubled[i,2,1] <> 0) AND (doubled [i,2,2] = 0) THEN BEGIN
    FOR k:= 1 to numind[i] DO
     IF pro[i,k] = 1 THEN  BEGIN
      pro[i,k]:=2;
      doubled[i,2,2] := k;
     END;
    END ELSE IF (doubled[i,2,2] <> 0) AND (doubled[i,2,1] = 0) THEN BEGIN
     FOR k:= 1 to numind[i] DO
      IF pro[i,k] = 1 THEN  BEGIN
       pro[i,k]:=2;
       doubled[i,2,1] := k;
      END;
    END;
 90:
 close(simped);
END; {readped}



PROCEDURE replicate;

VAR
 countp,a1,a2,number,i,j,k,ll,mm:integer;


 PROCEDURE founders; (* Assigns Genotypes to Founders *)
 VAR
  i,j,k,m:integer;
  pp,qq,st:real; {allele freq at previous, current locus}
  {st = factor s or t, depending on whether allele at previous locus was 1 or 2}
 BEGIN
  IF locustype[1]=1 THEN low:=2 ELSE low:=1;
  FOR i:=1 TO numpeds DO

  BEGIN
   FOR j:=1 TO numind[i] DO
    IF (pa[i,j]=0) and (pro[i,j] <= 1) THEN BEGIN
     simmed[i,j]:=true;
     FOR m:=1 TO 2 DO  {for m-th allele in individual, m = 1, 2}
      FOR k:=low TO numloci DO BEGIN
       if k=low then begin {first locus}
        if doub<genefreq[k,1] {allele 1 at k-th locus}
        then locus[i,j,k,m]:=1 else locus[i,j,k,m]:=2;
        {ped i, indiv j, locus k, allele m in an individual}
       end else begin {subsequent loci}
        pp:=genefreq[k-1,1];
        qq:=genefreq[k,1];
        st:=sqrt(corr2*pp*(1.0-pp)*qq*(1.0-qq));
        if locus[i,j,k-1,m]=1
         then if doub<(1.0-qq-st/pp)
          then locus[i,j,k,m]:=2
          else locus[i,j,k,m]:=1
         else if doub<(qq-st/(1.0-pp))
          then locus[i,j,k,m]:=1
          else locus[i,j,k,m]:=2;
       end; {else begin}
      END; {for k}
     IF (sexlink=1) AND (sex[i,j]=1) THEN
      FOR k:=low TO numloci DO  {changed 1 to low}
       locus[i,j,k,1]:=locus[i,j,k,2];
    END; {if}
  End; {if}

 END; {founders}


 PROCEDURE finishsimulation;
 VAR
  i,j,k,m,unsimmed,chro,ori:integer;
 BEGIN
  unsimmed:=999;
  FOR i:=1 TO numpeds DO
   simmed[i,0]:=true;
  WHILE unsimmed>0 DO BEGIN
   unsimmed:=0;
   FOR i:=1 TO numpeds DO
    FOR j:=1 TO numind[i] DO
     IF (pa[i,j] > 0) THEN BEGIN
      IF NOT simmed[i,j] THEN
       IF (NOT simmed[i,pa[i,j]]) OR (NOT simmed[i,ma[i,j]]) THEN
       unsimmed:=unsimmed+1 ELSE BEGIN
	FOR m:=1 TO 2 DO BEGIN
	 CASE m OF
	  1:ori:=pa[i,j];
	  2:ori:=ma[i,j];
	 END;
	 IF doub<0.5 THEN chro:=1 ELSE chro:=2;
	 FOR k:=low TO numloci DO BEGIN
	  locus[i,j,k,m]:=locus[i,ori,k,chro];
	  IF k<numloci THEN IF doub<theta[k,m]
	  THEN IF chro=1 THEN chro:=2 ELSE chro:=1
	 END;
	END;
	IF (sexlink=1) AND (sex[i,j]=1) THEN
	 FOR k:=low TO numloci DO locus[i,j,k,1]:=locus[i,j,k,2];
	simmed[i,j]:=true;
       END;
     END ELSE IF (pro[i,j]>1) THEN IF simmed[i,doubled[i,pro[i,j],2]] THEN BEGIN
       simmed[i,j]:=true;
       FOR k:= low TO numloci DO
	   FOR m:= 1 to 2 DO
	      locus[i,j,k,m]:=locus[i,doubled[i,pro[i,j],2],k,m];
     END ELSE IF NOT simmed[i,j] then unsimmed:=unsimmed+1;
  END;
 END; {finishsimulation}

BEGIN {replicate}
 founders;
 finishsimulation;
 FOR i:=1 TO numpeds DO
  FOR j:=1 TO numind[i] DO BEGIN
   number:=((ss-1)*numpeds)+i;
   write(pedfile,number:5,indnum[i,j]:4,pa[i,j]:4,ma[i,j]:4,fo[i,j]:4,
     nps[i,j]:4,nms[i,j]:4,sex[i,j]:3,pro[i,j]:3);
   IF locustype[1]=1 THEN BEGIN write(pedfile,disease[i,j]:6);
    IF liabnum[1]>1 THEN write(pedfile,liabcl[i,j]:3);
   END;
   FOR k:=low TO numloci DO
    IF typed[i,j,k]=1 THEN
     CASE locustype[k] OF
      0:BEGIN
	 IF locus[i,j,k,1]<=locus[i,j,k,2]
	 THEN BEGIN a1:=locus[i,j,k,1];a2:=locus[i,j,k,2] END
	 ELSE BEGIN a1:=locus[i,j,k,2];a2:=locus[i,j,k,1];END;
	 countp:=0;
	 FOR ll:=1 TO numall[k] DO
	  FOR mm:=ll TO numall[k] DO BEGIN
	   countp:=countp+1;
	   IF (a1=ll) AND (a2=mm) THEN BEGIN
	    IF (a1=a2)
	    THEN qphen[i,j,k]:=nordev(mean[k,countp],vari[k])
	    ELSE qphen[i,j,k]:=nordev(mean[k,countp],(vari[k]*multvar[k]));
	    IF pro[i,j] > 1 THEN IF (doubled[i,pro[i,j],1] < j) THEN
	       qphen[i,j,k]:=qphen[i,doubled[i,pro[i,j],1],k] ELSE
	       IF (doubled[i,pro[i,j],2] < j) THEN
	       qphen[i,j,k]:=qphen[i,doubled[i,pro[i,j],2],k];
	    write(pedfile,qphen[i,j,k]:12:6);
	   END;
	  END;
	END;
      1:BEGIN
	 IF locus[i,j,k,1]<=locus[i,j,k,2]
	 THEN BEGIN a1:=locus[i,j,k,1];a2:=locus[i,j,k,2];END
	 ELSE BEGIN a1:=locus[i,j,k,2];a2:=locus[i,j,k,1];END;
	 countp:=0;
	 FOR ll:=1 TO numall[k] DO
	  FOR mm:=ll TO numall[k] DO BEGIN
	   countp:=countp+1;
	   IF (a1=ll) AND (a2=mm) THEN BEGIN
	    IF doub<pen[k,countp] THEN aphen[i,j,k]:=2 ELSE aphen[i,j,k]:=1;
	    IF pro[i,j] > 1 THEN IF (doubled[i,pro[i,j],1] < j) THEN
	       aphen[i,j,k]:=aphen[i,doubled[i,pro[i,j],1],k] ELSE
	       IF (doubled[i,pro[i,j],2] < j) THEN
	       aphen[i,j,k]:=aphen[i,doubled[i,pro[i,j],2],k];
	    write(pedfile,aphen[i,j,k]:8);
	   END;
	  END;
	END;
      2:BEGIN
	 a1:=locus[i,j,k,1];
	 a2:=locus[i,j,k,2];
	 FOR ll:=1 TO nfact[k] DO
	  facttot[k,ll]:=binall[k,a1,ll];
	 FOR ll:=1 TO nfact[k] DO
	  IF binall[k,a2,ll]=1 THEN facttot[k,ll]:=binall[k,a2,ll];
	 write(pedfile,'   ');
	 FOR ll:=1 TO nfact[k] DO
	  write(pedfile,facttot[k,ll]:2);
	END;
      3:write(pedfile,locus[i,j,k,1]:5,locus[i,j,k,2]:3);
     END
    ELSE CASE locustype[k] OF 0:write(pedfile,'     0.0    ');
     1:write(pedfile,'       0');
     2:BEGIN
	write(pedfile,'   ');
	FOR ll:=1 TO nfact[k] DO write(pedfile,' 0');
       END;
     3:write(pedfile,'    0  0');
    END;
   writeln(pedfile);
   simmed[i,j]:=false;
  END;
END; {replicate}


PROCEDURE writesim;

VAR
 simout:text;
 jjj,cnt:integer;

BEGIN
 if verbose then writeln('Opening "simout.dat" file for output');
 assign(simout,'simout.dat');
 rewrite(simout);
 writeln(simout,' The random number seeds are: (obsolete)');
 writeln(simout,' The number of replicates is:',numreps:13);
 writeln(simout,' The requested proportion of unlinked families is:  0.000');
 writeln(simout,' The trait locus is number:   1');
 writeln(simout,'    Summary Statistics about simped.dat');
 writeln(simout,' Number of Pedigrees',numpeds:10);
 write(simout,' Number of People');
 cnt:=0;
 FOR jjj:=1 TO numpeds DO
  cnt:=cnt+numind[jjj];
 writeln(simout,cnt:13);
 writeln(simout,'LD (r^2) between consecutive SNPs =',corr2:7:4);
 close(simout);
END; {writesim}


procedure initialize;
var i1,i2,i3,i4:integer;
begin
 {writeln;}
 writeln('Program SIMULATE3 version ',version);
 if verbose then begin
  writeln;
  writeln('Input files required:');
  writeln('  PROBLEM.DAT  number of replicates');
  writeln('  SIMDATA.DAT  datafile in LINKAGE format');
  writeln('  SIMPED.DAT   pedfile in LINKAGE format');
  writeln;
  writeln('Maximum values of parameters:');
  writeln(maxpeds:8,' pedigrees');
  writeln(maxloci:8,' loci');
  writeln(maxind:8,' individuals in any one pedigree');
  writeln(maxall:8,' alleles at any one locus');
  writeln(maxfact:8,' binary factors');
  writeln;
 end;

 bad:=false;    {Initialize arrays}
 corr2:=0.0;
 for i1:=1 to maxpeds do for i2:=1 to 5 do for i3:=1 to 2 do
  doubled[i1,i2,i3]:=0;
 for i1:=1 to maxloci do begin
  liabnum[i1]:=0;
  locustype[i1]:=0;
  numall[i1]:=0;
 end;
 for i1:=0 to maxloci do for i2:=1 to 2 do
  theta[i1,i2]:=0;
 for i1:=1 to maxloci do for i2:=1 to maxall do
  genefreq[i1,i2]:=0;
 for i1:=1 to maxpeds do numind[i1]:=0;
 for i1:=1 to maxpeds do for i2:=1 to maxind do begin
  indnum[i1,i2]:=0;
  pa[i1,i2]:=0;
  ma[i1,i2]:=0;
  fo[i1,i2]:=0;
  nps[i1,i2]:=0;
  nms[i1,i2]:=0;
  pro[i1,i2]:=0;
  disease[i1,i2]:=0;
  liabcl[i1,i2]:=0;
 end;
 for i1:=1 to maxpeds do for i2:=1 to maxind do
  sex[i1,i2]:=1;
 for i1:=1 to maxpeds do for i2:=1 to maxind do for i3:=1 to maxloci do typed[i1,i2,i3]:=0;
 for i1:=1 to maxpeds do for i2:=1 to maxind do for i3:=1 to maxloci do
  for i4:=1 to 2 do locus[i1,i2,i3,i4]:=1;
 for i1:=1 to maxpeds do for i2:=0 to maxind do simmed[i1,i2]:=false;
 for i1:=1 to maxloci do begin
  nfact[i1]:=0;
  ntrait[i1]:=0;
  multvar[i1]:=0;
  vari[i1]:=0;
 end;
 for i1:=1 to maxloci do for i2:=1 to maxgen do begin
  mean[i1,i2]:=0;
  pen[i1,i2]:=0;
 end;
 for i1:=1 to maxloci do for i2:=1 to maxall do for i3:=1 to maxfact do binall[i1,i2,i3]:=0;
 for i1:=1 to maxloci do for i2:=1 to maxfact do facttot[i1,i2]:=0;
 for i1:=1 to maxpeds do for i2:=1 to maxind do for i3:=1 to maxloci do qphen[i1,i2,i3]:=0;
 for i1:=1 to maxpeds do for i2:=1 to maxind do for i3:=1 to maxloci do aphen[i1,i2,i3]:=0;
 seqrun:=-1; {sequential number of run}

 raninit;

end; {initialize}


var exitsave:pointer;
{$I errtrap.pas}




begin {simulate3}

 exitsave:=exitproc;
 exitproc:=@errtrap;

 initialize;
 readfile; (* To Read Information From a File Named PROBLEM.DAT *)
 readdata;      (* To Read Parameter File, SIMDATA.DAT *)
 readped(bad);  (* To Read Pedigree File, SIMPED.DAT *)
 if bad then goto 99;
(*
 writeln;
 writeln('Loci and interlocus male ( / female) recombination fractions:');
 for ss:= 1 to numloci do begin
  write(ss:3); {changed from locustype[ss]}
  if ss < numloci then write(theta[ss,1]:7:3);
  if (sexdiff > 0)  AND (ss < numloci) then write('/',theta[ss,2]:7:3);
 end;
 writeln;
*)
 if verbose then writeln('Opening "pedfile.dat" file for output');
 assign(pedfile,'pedfile.dat');
 rewrite(pedfile);

 if verbose then writeln('Generating replicates');
 for ss:= 1 to numreps do begin  {To Simulate One Replicate.  Write to pedfile}
  replicate;
 end;
 writesim;      (* Write SIMOUT.DAT *)

 writefile;  {Write new PROBLEM.DAT file}
 99:
 close(pedfile);
 endseed;
end.
