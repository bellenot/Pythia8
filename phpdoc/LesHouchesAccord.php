<html>
<head>
<title>Les Houches Accord</title>
</head>
<body>

<script language=javascript type=text/javascript>
function stopRKey(evt) {
var evt = (evt) ? evt : ((event) ? event : null);
var node = (evt.target) ? evt.target :((evt.srcElement) ? evt.srcElement : null);
if ((evt.keyCode == 13) && (node.type=="text"))
{return false;}
}

document.onkeypress = stopRKey;
</script>
<?php
if($_POST['saved'] == 1) {
if($_POST['filepath'] != "files/") {
echo "<font color='red'>SETTINGS SAVED TO FILE</font><br/><br/>"; }
else {
echo "<font color='red'>NO FILE SELECTED YET.. PLEASE DO SO </font><a href='SaveSettings.php'>HERE</a><br/><br/>"; }
}
?>

<form method='post' action='LesHouchesAccord.php'>

<h2>Les Houches Accord</h2>

The Les Houches Accord (LHA) for user processes [<a href="Bibliography.php" target="page">Boo01</a>] is the 
standard way to input parton-level information from a 
matrix-elements-based generator into PYTHIA. The conventions for 
which information should be stored has been defined in a Fortran context, 
as two commonblocks. Here a C++ equivalent is defined, as two separate 
classes.

<p/>
The <code>LHAinit</code> and <code>LHAevnt</code> classes are base classes, 
containing reading and printout functions, plus a pure virtual function 
each. Derived classes have to provide these two virtual functions to do 
the actual work. The existing derived classes are for reading information 
from a Les Houches Event File (LHEF) or from the respective Fortran 
commonblock.

<p/>
Normally, pointers to objects of the derived classes should be handed
in with the <?php $filepath = $_GET["filepath"];
echo "<a href='ProgramFlow.php?filepath=".$filepath."' target='page'>";?>
<code>pythia.init( LHAinit*, LHAevnt*)</code></a> method. However, 
with the LHEF format a filename can replace the two pointers, see below. 

<h3>Initialization</h3>

The <code>LHAinit</code> class stores information equivalent to the 
<code>/HEPRUP/</code> commonblock, as required to initialize the event 
generation chain. The main difference is that the vector container 
now allows a flexible number of subprocesses to be defined. For the 
rest, names have been modified, since the 6-character-limit does not 
apply, and variables have been regrouped for clarity, but nothing 
fundamental is changed.

<p/>
The pure virtual function <code>set()</code> has to be implemented in the 
derived class, to set relevant information when called. It should
return <code>false</code> if it fails to set the info.

<p/>
Inside <code>set()</code>, such information can be set by the following 
methods:
<p/><code>method&nbsp; </code><strong> beamA( identity, energy, pdfGroup, pdfSet) &nbsp;</strong> <br/>
sets the properties of the first incoming beam (cf. the Fortran
<code>IDBMUP(1), EBMUP(1), PDFGUP(1), PDFSUP(1)</code>), and similarly 
a <code>beamB</code> method exists. The parton distribution information 
defaults to zero, meaning that internal sets are used.
  

<p/><code>method&nbsp; </code><strong> strategy( strategyIn) &nbsp;</strong> <br/>
sets the event weighting and cross section strategy. The default, 
provided in the class constructor, is 3.
<br/><code>argument</code><strong> strategyIn </strong>  : 
chosen strategy (cf. <code>IDWTUP</code>; see [<a href="Bibliography.php" target="page">Sjo06</a>] 
section 9.9.1 for extensive comments).
<br/><code>argumentoption </code><strong> 1</strong> :  events come with non-negative weight, given in units 
of pb, with an average that converges towards the cross section of the
process. PYTHIA is in charge of the event mixing, i.e. for each new
try decides which process should be generated, and then decides whether
is should be kept based on a comparison with <code>xMax</code>.
Accepted events therefore have unit weight.  
<br/><code>argumentoption </code><strong> -1</strong> :  as option 1, except that cross sections can now be
negative and events after unweighting have weight +-1. You can use 
<?php $filepath = $_GET["filepath"];
echo "<a href='EventInformation.php?filepath=".$filepath."' target='page'>";?><code>pythia.info.weight()</code></a> 
to find the weight of the current event. A correct event mixing requires 
that a process that can take both signs should be split in two, one limited 
to positive or zero and the other to negative or zero values, with 
<code>xMax</code> chosen appropriately for the two.  
<br/><code>argumentoption </code><strong> 2</strong> :  events come with non-negative weight, in unspecified 
units, but such that <code>xMax</code> can be used to unweight the events
to unit weight. Again PYTHIA is in charge of the event mixing.
The total cross section of a process is stored in 
<code>xSec</code>.  
<br/><code>argumentoption </code><strong> -2</strong> :  as option 2, except that cross sections can now be
negative and events after unweighting have weight +-1. As for option -1
processes with indeterminate sign should be split in two.  
<br/><code>argumentoption </code><strong> 3</strong> :  events come with unit weight, and are thus accepted 
as is. The total cross section of the process is stored in 
<code>xSec</code>.  
<br/><code>argumentoption </code><strong> -3</strong> :  as option 3, except that events now come with weight 
+-1. Unlike options -1 and -2 processes with indeterminate sign need not be 
split in two, unless you intend to mix with internal PYTHIA processes 
(see below).  
<br/><code>argumentoption </code><strong> 4</strong> :  events come with non-negative weight, given in units 
of pb, with an average that converges towards the cross section of the
process, like for option 1. No attempt is made to unweight the events, 
however, but all are generated in full, and retain their original weight.
For consistency with normal PYTHIA units, the weight stored in 
<code>pythia.info.weight()</code> has been converted to mb, however.
  
<br/><code>argumentoption </code><strong> -4</strong> :  as option 4, except that events now can come 
either with positive or negative weights.   
<br/><b>Note 1</b>: if several processes have already been mixed and 
stored in a common event file, either LHEF or some private format, it 
would be problematical to read back events in a different order. Since it 
is then not feasible to let PYTHIA pick the next process type, strategies 
+-1 and +-2 would not work. Instead strategy 3 would be the recommended
choice, or -3 if negative-weight events are required.
<br/><b>Note 2</b>: it is possible to switch on internally implemented 
processes and have PYTHIA mix these with LHA ones according to their relative 
cross sections for strategies +-1, +-2 and 3. It does not work for strategy 
-3 unless the positive and negative sectors of the cross sections are in 
separate subprocesses (as must always be the case for -1 and -2), since 
otherwise the overall mixture of PYTHIA and LHA processes will be off. 
Mixing is not possible for strategies +-4, since the weighting strategy 
is unspecified. 
  
  

<p/><code>method&nbsp; </code><strong> process( idProcess, xSec, xErr, xMax) &nbsp;</strong> <br/>
sets info on an allowed process (cf. <code>LPRUP, XSECUP, XERRUP, 
XMAXUP</code>). 
Each new call will append one more entry to the list of processes.
The choice of strategy determines which quantities are mandatory:
<code>xSec</code> for strategies +-2 and +-3, 
<code>xErr</code> never, and
<code>xMax</code> for strategies +-1 and +-2. 
  

<br/><b>Note</b>: PYTHIA does not make active use of the (optional)
<code>xErr</code> values, but calculates a statistical cross section 
error based on the spread of event-to-event weights. This should work 
fine for strategy options +-1, but not for the others. Specifically, 
for options +-2 and +-3 the weight spread may well vanish, and anyway 
is likely to be an underestimate of the true error. If the author of the 
LHA input information does provide error information you may use that -
this information is displayed at initialization. If not, then a relative
error decreasing like <i>1/sqrt(n_acc)</i>, where <i>n_acc</i>     
is the number of accepted events, should offer a reasonable estimate.

<p/>
Information is handed back by the following methods:
<p/><code>method&nbsp; </code><strong> idBeamA(), eBeamA(), pdfGroupBeamA(), pdfSetBeamA() &nbsp;</strong> <br/>
and similarly with <i>A -> B</i>, for the two beam particles.
  
<p/><code>method&nbsp; </code><strong> strategy() &nbsp;</strong> <br/>
for the strategy choice.
  
<p/><code>method&nbsp; </code><strong> size() &nbsp;</strong> <br/>
for the number of subprocesses.
  
<p/><code>method&nbsp; </code><strong> idProcess(i), xSec(i), xErr(i), xMax(i) &nbsp;</strong> <br/>
for process <code>i</code> in the range <code>0 &lt;= i &lt; 
size()</code>.   
  

<p/>
The information can also be printed using the <code>list()</code>
method, e.g. <code>LHAinitObject->list()</code>.
This is automatically done by the <code>pythia.init</code> call.

<h3>Event input</h3>

The <code>LHAevnt</code> class stores information equivalent to the 
<code>/HEPEUP/</code> commonblock, as required to hand in the next 
parton-level configuration for complete event generation. The main 
difference is that the vector container now allows a flexible number 
of partons to be defined. For the rest, names have been modified, 
since the 6-character-limit does not apply, and variables have been 
regrouped for clarity, but nothing fundamental is changed.

<p/>
The LHA standard is based on Fortran arrays beginning with
index 1, and mother information is defined accordingly. In order to 
be compatible with this convention, the zeroth line of the C++ particle
array is kept empty, so that index 1 also here corresponds to the first
particle. One small incompatibility is that the <code>size()</code> 
method returns the full size of the particle array, including the 
empty zeroth line, and thus is one larger than the true number of 
particles (<code>NUP</code>). 

<p/>
The pure virtual function <code>set(idProcess)</code> has to be 
implemented in the derived class, to set relevant information when 
called. For strategy options +-1 and +-2 the input 
<code>idProcess</code> value specifies which process that should be 
generated by <code>set(...)</code>, while <code>idProcess</code> is 
irrelevant for strategies +-3 and +-4. The <code>set(...)</code> 
function should return <code>false</code> if it fails to set the info, 
i.e. normally that the supply of events in a file is exhausted. 
If so, no event is generated, and <code>pythia.next()</code> returns 
<code>false</code>. You can then interrogate 
<?php $filepath = $_GET["filepath"];
echo "<a href='EventInformation.php?filepath=".$filepath."' target='page'>";?><code>pythia.info.atEndOfFile()</code></a> 
to confirm that indeed the failure is caused in the <code>set(...)</code> 
function, and decide to break out of the event generation loop.

<p/>
Inside a normal <code>set(...)</code> call, information can be set by 
the following methods:
<p/><code>method&nbsp; </code><strong> process( idProcess, weight, scale, alphaQED, alphaQCD) &nbsp;</strong> <br/>
tells which kind of process occured, with what weight, at what scale, 
and which <i>alpha_EM</i> and <i>alpha_strong</i> were used
(cf. <code>IDPRUP, XWTGUP, SCALUP, AQEDUP, AQCDUP</code>). This method 
also resets the size of the particle list, and adds the empty zeroth 
line, so it has to be called before the <code>particle</code> method below.
  
<p/><code>method&nbsp; </code><strong> particle( id, status, mother1, mother2, colourTag1, colourTag2, p_x, p_y, p_z, e, m, tau, spin) &nbsp;</strong> <br/>
gives the properties of the next particle handed in (cf. <code>IDUP, ISTUP, 
MOTHUP(1,..), MOTHUP(2,..), ICOLUP(1,..), ICOLUP(2,..),  PUP(J,..), 
VTIMUP, SPINUP</code>) .
  

<p/>
Information is handed back by the following methods:
<p/><code>method&nbsp; </code><strong> idProcess(), weight(), scale(), alphaQED(), alphaQCD() &nbsp;</strong> <br/>
Note that the weight stored in <code>pythia.info.weight()</code> as a rule
is not the same as the above <code>weight()</code>: the method here gives
the value before unweighting while the one in <code>info</code> gives
the one after unweighting and thus normally is 1 or -1. Only with strategy
options +-4 would the value in <code>info</code> be the same as here, 
except for a conversion from pb to mb. 
  
<p/><code>method&nbsp; </code><strong> size() &nbsp;</strong> <br/>
for the size of the particle array, which is one larger than the number 
of particles in the event, since the zeroth entry is kept empty 
(see above).
  
<p/><code>method&nbsp; </code><strong> id(i), status(i), mother1(i), mother2(i), col1(i), col2(i),px(i), py(i), pz(i), e(i), m(i), tau(i), spin(i) &nbsp;</strong> <br/>
for particle <code>i</code> in the range 
<code>0 &lt;= i &lt; size()</code>. (But again note that 
<code>i = 0</code> is an empty line, so the true range begins at 1.)   
  

<p/>
In the LHEF description [<a href="Bibliography.php" target="page">Alw06</a>] an extension to 
include information on the parton densities of the colliding partons
is suggested. This optional further information can be set by
<p/><code>method&nbsp; </code><strong> pdf( id1, id2, x1, x2, scalePDF, xpdf1, xpdf2) &nbsp;</strong> <br/>
which gives the flavours , the <i>x</i> and the <ie>Q</i> scale 
(in GeV) at which the parton densities <i>x*f_i(x, Q)</i> have been
evaluated.
  

<p/>
This information is returned by the methods
<p/><code>method&nbsp; </code><strong> pdfIsSet(), id1(), id2(), x1(), x2(), scalePDF(), xpdf1(), xpdf2() &nbsp;</strong> <br/>
where the first one tells whether this optional information has been set
for the current event. (<code>pdf(...)</code> must be called after the
<code>process(...)</code> call of the event for this to work.)
  

<p/>
The information can also be printed using the <code>list()</code>
method, e.g. <code>LHAevntObject->list()</code>.
In cases where the <code>LHAevntObject</code> is not available to the
user, the <code>pythia.LHAevntList()</code> method can be used, which 
is a wrapper for the above. 

<p/>
The LHA expects the decay of resonances to be included as part of the 
hard process, i.e. if unstable particles are produced in a process then
their decays are also described. This includes <i>Z^0, W^+-, H^0</i> 
and other short-lived particles in models beyond the Standard Model.
Should this not be the case then PYTHIA will perform the decays of all
resonances it knows how to do, in the same way as for internal processes.
Note that you will be on slippery ground if you then restrict the decay of
these resonances to specific allowed channels since, if this is not what 
was intended, you will obtain the wrong cross section and potentially the
wrong mix of different event types. (Since the original intention is
unknown, the cross section will not be corrected for the fraction of
open channels, i.e. the procedure used for internal processes is not
applied in this case.) 

<h3>An interface to Les Houches Event Files</h3>

The LHEF standard [<a href="Bibliography.php" target="page">Alw06</a>] specifies a format where a single file 
packs initialization and event information. This has become the most 
frequently used procedure to process external parton-level events in Pythia. 
Therefore a special 
<?php $filepath = $_GET["filepath"];
echo "<a href='ProgramFlow.php?filepath=".$filepath."' target='page'>";?><code>pythia.init(fileName)</code></a>
initialization option exists, where the LHEF name is provided as single
input. Internally this name is then used to create instances of two derived 
classes, <code>LHAinitLHEF</code> and <code>LHAevntLHEF</code>. Both of
them are allowed to read from the same LHEF, first the former and then 
the latter.

<p/>
An example how to generate events from a LHEF is found in 
<code>main12.cc</code>. Note the use of 
<code>pythia.info.atEndOfFile()</code> to find out when the whole
LHEF has been processed. 

<p/>
To allow the sequential use of several event files the 
<code>init(...)</code> method has an optional second argument: 
<code>pythia.init(fileName, bool skipInit = false)</code>.
If called with this argument <code>true</code> then there will
be no initialization, except that the existing <code>LHAinitLHEF</code> 
and <code>LHAevntLHEF</code> classes will be deleted and replaced 
by ones pointing to the new file. It is assumed (but never checked) that 
the initialization information is identical, and that the new file simply
contains further events of exactly the same kind as the previous one.
An example of this possibility, and the option to mix with internal 
processes, is found in <code>main13.cc</code>.

<h3>A runtime Fortran interface</h3>

The runtime Fortran interface requires linking to an external Fortran
code. In order to avoid problems with unresolved external references
when this interface is not used, the code has been put in a separate
<code>LHAFortran.h</code> file, that is not included in any of the
other library files. Instead it should be included in the 
user-supplied main program, together with the implementation of two
methods below that call the Fortran program to do its part of the job. 

<p/>
The <code>LHAinitFortran</code> class derives from <code>LHAinit</code>. 
It reads initialization information from the LHA standard 
Fortran commonblock, assuming this commonblock behaves like an 
<code>extern "C" struct</code> named <code>heprup_</code>. (Note the final
underscore, to match how the gcc compiler internally names Fortran
files.) 

<p/>
Initialization is with
<pre>
    LHAinitFortran lhaInit();
</pre>
i.e. does not require any arguments. 

<p/>
The user has to supply an implementation of the <code>fillHepRup()</code>
method, that is to do the actual calling of the external Fortran routine(s)
that fills the <code>HEPRUP</code> commonblock. The translation of this
information to the C++ structure is provided by the existing 
<code>set()</code> code. 

<p/>
The <code>LHAevntFortran</code> class derives from <code>LHAevnt</code>. 
It reads information on the next event, stored in the LHA 
standard Fortran commonblock, assuming this commonblock behaves like 
an <code>extern "C" struct</code> named <code>hepeup_</code>.

<p/>
Initialization is with
<pre>
    LHAevntFortran lhaEvnt();
</pre>
i.e. does not require any arguments. 

<p/>
The user has to supply an implementation of the <code>fillHepEup()</code>
method, that is to do the actual calling of the external Fortran routine(s)
that fills the <code>HEPEUP</code> commonblock. The translation of this
information to the C++ structure is provided by the existing 
<code>set()</code> code. 

<p/>
See further 
<?php $filepath = $_GET["filepath"];
echo "<a href='AccessPYTHIA6Processes.php?filepath=".$filepath."' target='page'>";?>here</a> for 
information how PYTHIA 6.4 can be linked to make use of this facility. 

</body>
</html>

<!-- Copyright (C) 2007 Torbjorn Sjostrand -->
