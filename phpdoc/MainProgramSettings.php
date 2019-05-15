<html>
<head>
<title>Main-Program Settings</title>
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

<form method='post' action='MainProgramSettings.php'>

<h2>Main-Program Settings</h2>

<h3>Introduction</h3>

The main program is up to the user to write. However, sample main
programs are provided, as documented further down on this page. In one 
such class of programs, key settings of the run are read in from a 
"cards file". These commands may be of two types<br/>
(a) instructions directly to <code>Pythia</code>, like which 
processes to generate, and<br/>
(b) instructions to the main program for what it should do, 
like how many events to generate, i.e. how many times 
<code>pythia.next()</code> should be called.<br/>
In principle these two kinds could be kept completely separate. 
However, to make life simpler, a number of useful main-program 
settings are defined on this page, so that they are recognized by 
the <code>Settings</code> machinery. They can thus be put among 
the other cards without distinction. It is up to you to decide which 
ones, if any, you actually want to use when you write your main program.
For convenience, the ones in the first section below, and some in the
second section, can also be interpreted directly by <code>Pythia</code>, 
while the subsequent ones really have to be used in your main program. 

<p/>
Once you have used the <code>pythia.readFile(fileName)</code> method to
read in the cards file, you can interrogate the <code>Settings</code>
database to make the values available in your main program. A slight
complication is that you need to use a different  <code>Settings</code>
method for each of the four possible return types that you want to 
extract. To save some typing the same method names are found directly 
in the <code>Pythia</code> class, and just send on to the
<code>Settings</code> ones to do the job, e.g.
<pre>
  bool   showCS = pythia.flag("Main:showChangedSettings");
  int    nEvent = pythia.mode("Main:numberOfEvents");
  double eCM    = pythia.parm("Main:eCM");
  string file   = pythia.word("Main:allSettingsFile"); 
</pre>

<h3>Incoming beams</h3>

Normally the identities and energies of the two incoming beam particles 
are given by the arguments of the various forms of the <code>init</code> 
call. These settings can be stored in an input "cards" file, in the 
following variables, and thereafter read in the user-written main program. 
As a shortcut, an <code>init()</code> call with no arguments will make use 
of the beam values directly. That is, if nothing is set, you will default 
to LHC at the nominal energy. Usage is purely optional. 

<br/><br/><table><tr><td><strong>Main:idBeamA  </td><td></td><td> <input type="text" name="1" value="2212" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2212</strong></code>)</td></tr></table>
The PDG <code>id</code> code for the first incoming particle.
</modeopen>

<br/><br/><table><tr><td><strong>Main:idBeamB  </td><td></td><td> <input type="text" name="2" value="2212" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2212</strong></code>)</td></tr></table>
The PDG <code>id</code> code for the second incoming particle.
</modeopen>

<br/><br/><strong>Main:inCMframe</strong>  <input type="radio" name="3" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="3" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Assume collisions occur in the CM frame.
  

<br/><br/><table><tr><td><strong>Main:eCM </td><td></td><td> <input type="text" name="4" value="14000." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>14000.</strong></code>; <code>minimum = 10.</code>)</td></tr></table>
Collision CM energy, to be given if <code>Main:inCMframe</code> is on. 
  

<br/><br/><table><tr><td><strong>Main:eBeamA </td><td></td><td> <input type="text" name="5" value="7000." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>7000.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The energy of the first incoming particle, moving in the 
<i>+z </i>direction, to be given if <code>Main:inCMframe</code> 
is off. If the particle energy is smaller than its mass
it is assumed to be at rest. 
  

<br/><br/><table><tr><td><strong>Main:eBeamB </td><td></td><td> <input type="text" name="6" value="7000." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>7000.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The energy of the second incoming particle, moving in the 
<i>-z</i> direction, to be given if <code>Main:inCMframe</code> 
is off. If the particle energy is smaller than its mass
it is assumed to be at rest.
  

<br/><br/><table><tr><td><strong>Main:LHEF  </td><td></td><td> <input type="text" name="7" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
The name of a Les Houches Event File. If you initialize with 
<code>init()</code> without any arguments, and <code>Main:LHEF</code>
has been set differently from its default value <code>void</code>, the 
initialization and subsequent run is based on the information stored 
in this file, overriding the beam-parameter input above.
  

<p/>
Currently there are no provisions for arbitrary beam directions,
but you can always rotate and boost the final 
<?php $filepath = $_GET["filepath"];
echo "<a href='EventRecord.php?filepath=".$filepath."' target='page'>";?>event record</a> appropriately.
For instance, consider two beams of equal energy but with a slight 
acollinearity: they are both an angle <i>chi</i> away from 
the <i>+-z</i> axis in the  <i>+x</i> direction, such that the
total acollinearity is <i>2 chi</i>. Then a boost <i>beta_x = chi</i>, 
achieved by <code>pythia.event.bst( chi, 0., 0.)</code>, moves
the event to the correct frame. 

<h3>Subruns</h3>

You can use <?php $filepath = $_GET["filepath"];
echo "<a href='ProgramFlow.php?filepath=".$filepath."' target='page'>";?>subruns</a> to carry out
several tasks in the same run. In that case you will need repeated
instances of the first setting below in your command file, and could
additionally use the second and third as well.

<br/><br/><table><tr><td><strong>Main:subrun  </td><td></td><td> <input type="text" name="8" value="-999" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-999</strong></code>; <code>minimum = 0</code>)</td></tr></table>
The number of the current subrun, a non-negative integer, put as
first line in a section of lines to be read for this particular subrun.
</modeopen>

<br/><br/><strong>Main:LHEFskipInit</strong>  <input type="radio" name="9" value="on"><strong>On</strong>
<input type="radio" name="9" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
If you read several Les Houches Event Files that you want to see 
considered as one single combined event sample you can set this flag
<code>on</code> after the first subrun to skip (most of) the  
(re-)initialization step.
  

<br/><br/><table><tr><td><strong>Main:numberOfSubruns  </td><td></td><td> <input type="text" name="10" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>)</td></tr></table>
The number of subruns you intend to use in your current run.  
Unlike the two settings above, <code>Pythia</code> itself will not
intepret this number, but you could e.g. have a loop in your main
program to loop over subruns from 0 through 
<code>numberOfSubruns - 1</code>. 
  

<h3>Run settings</h3>

Here further settings related to how many events to generate and whether
to print some information on data used in run. Again these variables 
can be set in an input "cards" file, and thereafter read out an used 
in the user-written main program. Usage is purely optional, but may help
you reduce the need to recompile your main program. 

<br/><br/><table><tr><td><strong>Main:numberOfEvents  </td><td></td><td> <input type="text" name="11" value="1000" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1000</strong></code>; <code>minimum = 0</code>)</td></tr></table>
The number of events to be generated.
</modeopen>

<br/><br/><table><tr><td><strong>Main:numberToList  </td><td></td><td> <input type="text" name="12" value="2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2</strong></code>; <code>minimum = 0</code>)</td></tr></table>
The number of events to list.
</modeopen>

<br/><br/><table><tr><td><strong>Main:timesToShow  </td><td></td><td> <input type="text" name="13" value="50" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>50</strong></code>; <code>minimum = 0</code>)</td></tr></table>
Print the number of events generated so far, this many times, 
i.e. once every <code>numberOfEvents/numberToShow</code> events.
</modeopen>

<br/><br/><table><tr><td><strong>Main:timesAllowErrors  </td><td></td><td> <input type="text" name="14" value="10" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>10</strong></code>)</td></tr></table>
Allow this many times that <code>pythia.next()</code> returns false, 
i.e. that an event is flawed, before aborting the run.
</modeopen>

<br/><br/><strong>Main:showChangedSettings</strong>  <input type="radio" name="15" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="15" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Print a list of the changed flag/mode/parameter/word settings.
  

<br/><br/><strong>Main:showAllSettings</strong>  <input type="radio" name="16" value="on"><strong>On</strong>
<input type="radio" name="16" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Print a list of all flag/mode/parameter/word settings.
  

<br/><br/><strong>Main:showChangedParticleData</strong>  <input type="radio" name="17" value="on"><strong>On</strong>
<input type="radio" name="17" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Print a list of particle and decay data for those particles 
that were changed (one way or another).
  

<br/><br/><strong>Main:showAllParticleData</strong>  <input type="radio" name="18" value="on"><strong>On</strong>
<input type="radio" name="18" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Print a list of all particle and decay data.
  

<br/><br/><strong>Main:writeChangedSettings</strong>  <input type="radio" name="19" value="on"><strong>On</strong>
<input type="radio" name="19" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Write a file with the changed flag/mode/parameter/word settings, in
a format appropriate to be read in at the beginning of a new  
run, using the <code>pythia.readFile(fileName)</code> method. 
  

<br/><br/><table><tr><td><strong>Main:changedSettingsFile  </td><td></td><td> <input type="text" name="20" value="currentSettings.cmnd" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>currentSettings.cmnd</strong></code>)</td></tr></table>
The name of the file to which the changed flag/mode/parameter/word
settings are written if <code>Main:writeChangedSettings</code>
is on. 
  

<br/><br/><strong>Main:writeAllSettings</strong>  <input type="radio" name="21" value="on"><strong>On</strong>
<input type="radio" name="21" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Write a file with all flag/mode/parameter/word settings, in
a format appropriate to be read in at the beginning of a new  
run, using the <code>pythia.readFile(fileName)</code> method. 
  

<br/><br/><table><tr><td><strong>Main:allSettingsFile  </td><td></td><td> <input type="text" name="22" value="allSettings.cmnd" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>allSettings.cmnd</strong></code>)</td></tr></table>
The name of the file to which a flag/mode/parameter/word 
settings are written if <code>Main:writeAllSettings</code>
is on. 
  

<br/><br/><strong>Main:showAllStatistics</strong>  <input type="radio" name="23" value="on"><strong>On</strong>
<input type="radio" name="23" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Print all available statistics or only the minimal set at the end 
of the run.
  

<h3>Sample main programs</h3>

To help exemplify what a main program could look like, a few simple
examples are provided: 

<ul>

<li><code>main01.cc</code> : a simple study of the charged multiplicity
for jet events at the LHC. (Brief example given in talks.)</li>

<li><code>main02.cc</code> : a simple study of the <i>pT</i> spectrum 
of Z bosons at the Tevatron. (Brief example given in talks.)</li>

<li><code>main03.cc</code> : a simple single-particle analysis of jet 
events, where input is set by <code>main03.cmnd</code> "cards file".</li>

<li><code>main04.cc</code> : a simple study of several different kinds 
of events, with the choice to be made in the <code>main04.cmnd</code> 
"cards file".</li>

<li><code>main05.cc</code> : generation of QCD jet events at the LHC, 
with jet analysis using the <code>CellJet</code> cone-jet finder.</li>

<li><code>main06.cc</code> : tests of cross sections for elastic and 
diffractive topologies, using <code>main06.cmnd</code> to pick process.</li>

<li><code>main07.cc</code> : tests of cross sections for minimum-bias 
events, using <code>main07.cmnd</code> to pick options.</li>

<li><code>main08.cc</code> : generation of the QCD jet cross section
by splitting the run into subruns, each in its own <i>pT</i> bin,
and adding the results properly reweighted. Two options, with limits 
set either in the main program or by subrun specification in the
<code>main08.cmnd</code> file.</li>

<li><code>main09.cc</code> : generation of LEP1 hadronic events, i.e. 
<i>e^+e^- -> gamma*/Z^0 -> q qbar</i>, with charged multiplicity, 
sphericity, thrust and jet analysis.</li>

<li><code>main10.cc</code> : illustration how userHooks can be used
interact directly with the event-generation process.</li>

<li><code>main11.cc</code> : generation of two predetermined hard
interactions in each event.</li>

<li><code>main12.cc</code> : a study of top events, fed in from the 
Les Houches Event File <code>ttbar.lhe</code>, here generated by 
<code>main53.f</code>. This file currently only contains 100 events 
so as not to make the distributed PYTHIA package too big, and so serves 
mainly as a demonstration of the principles involved. </li> 

<li><code>main13.cc</code> : a systematic comparison of several 
cross section values with their corresponding values in PYTHIA 6.4,
the latter available as a table in the code.

<li><code>main21.cc</code> : an example how parton-level configurations
can be input directly for hadronization, without being tied to the
full process-generation machinery, e.g. to study the hadronization of 
junction topologies.</li>

<li><code>main22.cc</code> : tests of internally implemented cross sections
for Supersymmetric particle production, with SYSY spectrum defined in
<code>main22.spc</code> and settings in <code>main22.cmnd</code>.</li>

<li><code>main23.cc</code> : shows how an external decay handler can 
be linked to handle the decays of some particles.</li>

<li><code>main24.cc</code> : shows how an external random number 
generator can be linked to replace the internal one.</li>

<li><code>main25.cc</code> : shows how an external process can be 
implemented as a new class derived from a PYTHIA base class, and then
handed in for generation as with a normal internal process.</li>

<li><code>main31.cc</code> : similar to main01, except that the 
event record is output in the HepMC event record format. Requires 
that HepMC and CLHEP are properly linked.</li>

<li><code>main32.cc</code> : a streamlined version, where all input on
the event generation is taken from the <code>main32.cmnd</code> file,
and the ony thing done is to generate events and store them in a
file in the HepMC event record format. All physics studies will have
to be done afterwards. Requires that HepMC and CLHEP are properly 
linked.</li>

<li><code>main41.cc</code> : a test of the shape of parton densities,
as a check prior to using a given PDF set in a generator.  Requires 
that LHAPDF is properly linked.</li>

<li><code>main42.cc</code> : compares the charged multiplicity 
distribution, and a few other aspects, between default PYTHIA PDF and 
another one. Requires that LHAPDF is properly linked.</li>

<li><code>main51.cc</code> : a simple example how the Les Houches
Accord interface, plus a few more Fortran-to-C++ commands, allows
hard processes to be generated by PYTHIA 6.4 and then processed 
further by PYTHIA 8. Requires that PYTHIA 6.4 is properly linked.</li>

<li><code>main52.cc</code> : a fairly extensive study of 
event properties, with hard processes generated by PYTHIA 6.4. 
It reads in a <code>main52.fcmnd</code> file with commands specfically
for the Fortran PYTHIA 6.4 program and another <code>main52.ccmnd</code> 
file illustrating several of the settings listed on these pages.
Requires that PYTHIA 6.4 is properly linked.</li>

<li><code>main53.f</code> : a Fortran program (!) showing how 
PYTHIA 6.4 can be used to generate a Les Houches Event File 
<code>ttbar.lhe</code> with top events (which is used as input by
<code>main12.cc</code>). This program can easily be modified to 
generate other files, bigger and/or for other processes.
Requires that PYTHIA 6.4 is properly linked.</li>

<li><code>main54.cc</code> : a final example where PYTHIA 6.4 is used 
to generate hard processes, which are directly input to be generated
in full by the internal machinery, using the settings in 
<code>main54.cmnd</code>, and the output consists of a file with 
HepMC event records for further analysis. Requires that PYTHIA 6.4, 
HepMC and CLHEP are properly linked.</li>

</ul>

<h3>Spares</h3>

For currently unforeseen purposes, a few dummy settings are made 
available here. The user can set the desired value in a "cards file"
and then use that value in the main program as desired.

<br/><br/><strong>Main:spareFlag1</strong>  <input type="radio" name="24" value="on"><strong>On</strong>
<input type="radio" name="24" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
  

<br/><br/><strong>Main:spareFlag2</strong>  <input type="radio" name="25" value="on"><strong>On</strong>
<input type="radio" name="25" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
  

<br/><br/><strong>Main:spareFlag3</strong>  <input type="radio" name="26" value="on"><strong>On</strong>
<input type="radio" name="26" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
  

<br/><br/><table><tr><td><strong>Main:spareMode1  </td><td></td><td> <input type="text" name="27" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>)</td></tr></table>
</modeopen>

<br/><br/><table><tr><td><strong>Main:spareMode2  </td><td></td><td> <input type="text" name="28" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>)</td></tr></table>
</modeopen>

<br/><br/><table><tr><td><strong>Main:spareMode3  </td><td></td><td> <input type="text" name="29" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>)</td></tr></table>
</modeopen>

<br/><br/><table><tr><td><strong>Main:spareParm1 </td><td></td><td> <input type="text" name="30" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
  

<br/><br/><table><tr><td><strong>Main:spareParm2 </td><td></td><td> <input type="text" name="31" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
  

<br/><br/><table><tr><td><strong>Main:spareParm3 </td><td></td><td> <input type="text" name="32" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
  

<br/><br/><table><tr><td><strong>Main:spareWord1  </td><td></td><td> <input type="text" name="33" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
  

<br/><br/><table><tr><td><strong>Main:spareWord2  </td><td></td><td> <input type="text" name="34" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
  

<br/><br/><table><tr><td><strong>Main:spareWord3  </td><td></td><td> <input type="text" name="35" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
  

<input type="hidden" name="saved" value="1"/>

<?php
echo "<input type='hidden' name='filepath' value='".$_GET["filepath"]."'/>"?>

<table width="100%"><tr><td align="right"><input type="submit" value="Save Settings" /></td></tr></table>
</form>

<?php

if($_POST["saved"] == 1)
{
$filepath = $_POST["filepath"];
$handle = fopen($filepath, 'a');

if($_POST["1"] != "2212")
{
$data = "Main:idBeamA = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "2212")
{
$data = "Main:idBeamB = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "on")
{
$data = "Main:inCMframe = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "14000.")
{
$data = "Main:eCM = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "7000.")
{
$data = "Main:eBeamA = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "7000.")
{
$data = "Main:eBeamB = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "void")
{
$data = "Main:LHEF = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "-999")
{
$data = "Main:subrun = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "off")
{
$data = "Main:LHEFskipInit = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "0")
{
$data = "Main:numberOfSubruns = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "1000")
{
$data = "Main:numberOfEvents = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "2")
{
$data = "Main:numberToList = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "50")
{
$data = "Main:timesToShow = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "10")
{
$data = "Main:timesAllowErrors = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "on")
{
$data = "Main:showChangedSettings = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "off")
{
$data = "Main:showAllSettings = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "off")
{
$data = "Main:showChangedParticleData = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "off")
{
$data = "Main:showAllParticleData = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "off")
{
$data = "Main:writeChangedSettings = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "currentSettings.cmnd")
{
$data = "Main:changedSettingsFile = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "off")
{
$data = "Main:writeAllSettings = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "allSettings.cmnd")
{
$data = "Main:allSettingsFile = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "off")
{
$data = "Main:showAllStatistics = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "off")
{
$data = "Main:spareFlag1 = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "off")
{
$data = "Main:spareFlag2 = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
if($_POST["26"] != "off")
{
$data = "Main:spareFlag3 = ".$_POST["26"]."\n";
fwrite($handle,$data);
}
if($_POST["27"] != "0")
{
$data = "Main:spareMode1 = ".$_POST["27"]."\n";
fwrite($handle,$data);
}
if($_POST["28"] != "0")
{
$data = "Main:spareMode2 = ".$_POST["28"]."\n";
fwrite($handle,$data);
}
if($_POST["29"] != "0")
{
$data = "Main:spareMode3 = ".$_POST["29"]."\n";
fwrite($handle,$data);
}
if($_POST["30"] != "0.")
{
$data = "Main:spareParm1 = ".$_POST["30"]."\n";
fwrite($handle,$data);
}
if($_POST["31"] != "0.")
{
$data = "Main:spareParm2 = ".$_POST["31"]."\n";
fwrite($handle,$data);
}
if($_POST["32"] != "0.")
{
$data = "Main:spareParm3 = ".$_POST["32"]."\n";
fwrite($handle,$data);
}
if($_POST["33"] != "void")
{
$data = "Main:spareWord1 = ".$_POST["33"]."\n";
fwrite($handle,$data);
}
if($_POST["34"] != "void")
{
$data = "Main:spareWord2 = ".$_POST["34"]."\n";
fwrite($handle,$data);
}
if($_POST["35"] != "void")
{
$data = "Main:spareWord3 = ".$_POST["35"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2007 Torbjorn Sjostrand -->
