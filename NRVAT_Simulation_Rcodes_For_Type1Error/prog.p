{procedures and functions}


 FUNCTION iint64():QWord;
 {
  This function returns a 64 but integer random number when the
  values of u,v,w have already been set by the raninit procedure.
  Repeatedly calling this function will get many random numbers.

  Based on Press WH, Teukolsky SA, Vetterling WT, Flannery BP (2007)
  "Numerical recipes 3rd edition: The art of scientific computing."
  Cambridge University Press, Cambridge, UK; New York
 }
 CONST
  LL28: qword = 2862933555777941757;
  LL70: qword = 7046029254386353087;
  U42: longword = 4294957665;
 VAR x: QWord;
 BEGIN {iint64}
  u := u * LL28 + LL70;
  v := v xor (v shr 17); v := v xor (v shl 31); v := v xor (v shr 8);
  { $ffffffff = %11111111111111111111111111111111 is hexadecimal for 4294967295 }
  w := U42 * (w and $ffffffff) + (w shr 32);
  x := u xor (u shl 21);  x := x xor (x shr 35);  x := x xor (x shl 4);
  iint64 := (x + v) xor w;
 END; {iint64}



 FUNCTION doub():double;
 {
  This function returns a double precision floating random
  number in the range 0 to 1 when the values of u,v,w have
  already been set by the raninit procedure. Repeatedly
  calling this function will get many random numbers.
 }
 CONST rr: double = 5.42101086242752217E-20;
 BEGIN {doub}
  doub := rr * iint64;
 END; {doub}



 PROCEDURE raninit;
 {
  This function initiates the random number generator, iint64, that is,
  it sets the states of global variables u,v,w, given a seed j ("seed").
  Note: int64 is a predefined number type in Pascal; need to name desired
  function iint64 instead. Generator has period 3.138 x 10^57.

  Based on Press WH, Teukolsky SA, Vetterling WT, Flannery BP (2007)
  "Numerical recipes 3rd edition: The art of scientific computing."
  Cambridge University Press, Cambridge, UK; New York

  Required global declarations:
 VAR
  u,v,w:QWord;
  seed:QWord; Called j in Press book.
 }

 CONST LL41: qword = 4101842887655102017;

 VAR
  seedfile:text;
  ok:boolean;
  rr:double;
  res:longword;

 BEGIN {raninit}
  assign(seedfile,'seed.txt');
  {writeln('Opening "seed.txt" file to read positive integer seed');}
  reset(seedfile);
  read(seedfile,rr);
  if (rr<1.0) then begin
   writeln('ERROR: Random seed must be a positive integer number');
   halt;
  end;
  reset(seedfile);
  read(seedfile,seed);
  close(seedfile);
  {writeln('Initial seed = ',seed);}

  if (seed=LL41) then begin
   writeln('ERROR: The seed must not be equal to the following number:');
   writeln(' ',LL41);
   halt;
  end;

  {Check proper sizes of numbers}
  ok := (sizeof(res)=4) and (sizeof(seed)=8);
  if not ok then begin
   writeln('ERROR: One or both of the following number sizes are incorrect:');
   writeln(' longword = 4, qword = 8. Random number generator must be updated.');
   halt;
  end;

  v := LL41;
  w := 1;
  u := seed xor v; iint64;
  v := u; iint64;
  w := v; iint64;
 END; {raninit}



 procedure endseed;
 {write new seed file with last used random number}
 VAR
  seedfile:text;
  ss:qword;
 BEGIN {endseed}
  assign(seedfile,'seed.txt');
  rewrite(seedfile);
  ss:=iint64;
  writeln(seedfile,' ',ss);
  close(seedfile);
  {writeln('Updated seed written to "seed.txt" file: ',ss);}
 end; {endseed}






 procedure increase(var rrval:real);
 {Increase argument by 1.0}
 begin
  rrval:=rrval+1.0;
 end;


