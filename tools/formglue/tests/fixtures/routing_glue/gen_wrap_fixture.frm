#-
* Genuine #write wrap fixture: 78-char symbol names force
* BOTH wrap styles: backslash mid-word splits + bare indented splits.
* (102-char names crash FORM 5.0.1 #write with rc=134.)
Off statistics;
Symbols
  aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaax1,
  bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbx2,
  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccx3,
  ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddx4;
L FIX = (
   aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaax1
  +bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbx2
  +ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccx3
  +ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddx4)^6;
.sort
#$n = termsin_(FIX);
#write "COUNT fixture_terms = `$n'"
Format nospaces;
#write <wrap_fixture.out> "FIX=%E;", FIX
#close <wrap_fixture.out>
.end
