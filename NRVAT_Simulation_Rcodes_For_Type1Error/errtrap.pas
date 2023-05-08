procedure errtrap;

(*Error trapping routine for Borland Pascal 7.0.  See Language guide p. 301
  Requires outfile:text.  Prompts user to press <Enter>.

 Program will stop with the following exitcodes:
   error                     errorlevel
   ------------------------------------
   none                           0
   user interrupt                 1
   run time         exit code generated
   ------------------------------------

 Usage:

 In main program, declare variable exitsave:pointer.
 First two statements in main program:
  exitsave:=exitproc;
  exitproc:=@errtrap;

Thus, the main program would start as follows:

-------------
var exitsave:pointer;
{$I errtrap.pas}

begin {main}

 exitsave:=exitproc;
 exitproc:=@errtrap;
-------------
*)


Begin
 if erroraddr<>nil then begin
  write(' **ERROR: ');
  case exitcode of
     2:writeln('File not found');
     3:writeln('Path not found');
     4:writeln('Too many open files');
     5:writeln('File access denied');
     6:writeln('Invalid file handle');
    12:writeln('Invalid file access code');
    15:writeln('Invalid drive number');
    16:writeln('Cannot remove current directory');
    17:writeln('Cannot rename across drives');
    18:writeln('No more files.');
   100:writeln('Disk read error');
   101:writeln('Disk write error');
   102:writeln('File not assigned');
   103:writeln('File not open');
   104:writeln('File not open for input');
   105:writeln('File not open for output');
   106:writeln('Invalid numeric format');
   150:writeln('Disk is write protected');
   151:writeln('Unknown unit');
   152:writeln('Drive not ready');
   153:writeln('Unknown command');
   154:writeln('CRC error in data');
   155:writeln('Bad drive request structure length');
   156:writeln('Disk seek error');
   157:writeln('Unknown media type');
   158:writeln('Sector not found');
   159:writeln('Printer out of paper');
   160:writeln('Device write fault');
   161:writeln('Device read fault');
   162:writeln('Hardware failure');
   200:writeln('Division by zero');
   201:writeln('Range check error');
   202:writeln('Stack overflow.  Probably due to an undeclared loop.');
   203:writeln('Heap overflow.  Insufficient free memory.');
   204:writeln('Invalid pointer operation');
   205:writeln('Floating point overflow');
   206:writeln('Floating point underflow');
   207:writeln('Invalid floating point operation');
   208:writeln('Overlay manager not installed');
   209:writeln('Overlay file read error');
   211:writeln('Call to abstract method');
   212:writeln('Stream registration error');
   213:writeln('Collection index out of range');
   214:writeln('Collection overflow error');
   216:writeln('General protection fault (protected mode)');
  else writeln('Unspecified error');
  end;
  erroraddr:=nil;  {This suppresses error report by system}
 end;
 if exitcode=255 then begin
  exitcode:=1;
  writeln(' **User interrupt**');
 end;
 exitproc:=exitsave;
end; {errtrap}
