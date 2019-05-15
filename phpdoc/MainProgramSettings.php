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
programs are provided. In one such class of programs, key settings 
of the run are read in from a "cards file".
These commands may be of two types<br/>
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

<p/>
Once you have used the <code>pythia.readFile(...)</code> method to
read in the cards file, you can interrogate the <code>Settings</code>
database to make the values available in your main program. A slight
complication is that you need to use a different method for each of 
the four possible return types that you want to extract, e.g.:
<pre>
  bool   showCS = pythia.settings.flag("Main:showChangedSettings");
  int    nEvent = pythia.settings.mode("Main:numberOfEvents");
  double eCM    = pythia.settings.parm("Main:eCM");
  string file   = pythia.settings.word("Main:allSettingsFile"); 
</pre>
To save some typing, the same method names are found directly in the
<code>Pythia</code> class, and just send on to <code>Settings</code>
to do the job, so that you can use <code>pythia.flag(...)</code> 
instead of <code>pythia.settings.flag(...)</code>, etc.

<h3>Incoming beams</h3>

Normally the identities and energies of the two incoming beam particles 
are given by the arguments of the init call. These settings can be
stored in an input "cards" file, in the following variables, and 
thereafter read in the user-written main program. Usage is purely 
optional. 

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
  

<br/><br/><table><tr><td><strong>Main:eCM </td><td></td><td> <input type="text" name="4" value="1960." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1960.</strong></code>; <code>minimum = 10.</code>)</td></tr></table>
Collision CM energy, to be given if <code>Main:inCMframe</code> is on. 
  

<br/><br/><table><tr><td><strong>Main:eBeamA </td><td></td><td> <input type="text" name="5" value="7000." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>7000.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The energy of the first incoming particle, moving in the 
<i>+z </i>direction. If the particle energy is smaller than its mass
it is assumed to be at rest. 
  

<br/><br/><table><tr><td><strong>Main:eBeamB </td><td></td><td> <input type="text" name="6" value="7000." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>7000.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
The energy of the second incoming particle, moving in the 
<i>-z</i> direction. If the particle energy is smaller than its mass
it is assumed to be at rest.
  

<h3>Run settings</h3>

Here further settings related to how many events to generate and whether
to print some information on data used in run. Again these variables 
can be set in an input "cards" file, and thereafter read out an used 
in the user-written main program. Usage is purely optional. 

<br/><br/><table><tr><td><strong>Main:numberOfEvents  </td><td></td><td> <input type="text" name="7" value="1000" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1000</strong></code>; <code>minimum = 0</code>)</td></tr></table>
The number of events to be generated.
</modeopen>

<br/><br/><table><tr><td><strong>Main:numberToList  </td><td></td><td> <input type="text" name="8" value="2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2</strong></code>; <code>minimum = 0</code>)</td></tr></table>
The number of events to list.
</modeopen>

<br/><br/><table><tr><td><strong>Main:timesToShow  </td><td></td><td> <input type="text" name="9" value="50" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>50</strong></code>; <code>minimum = 0</code>)</td></tr></table>
Print the number of events generated so far, this many times, 
i.e. once every <code>numberOfEvents/numberToShow</code> events.
</modeopen>

<br/><br/><table><tr><td><strong>Main:timesAllowErrors  </td><td></td><td> <input type="text" name="10" value="10" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>10</strong></code>)</td></tr></table>
Allow this many times that <code>pythia.next()</code> returns false, 
i.e. that an event is flawed, before aborting the run.
</modeopen>

<br/><br/><strong>Main:showChangedSettings</strong>  <input type="radio" name="11" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="11" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Print a list of the changed flag/mode/parameter/word settings.
  

<br/><br/><strong>Main:showAllSettings</strong>  <input type="radio" name="12" value="on"><strong>On</strong>
<input type="radio" name="12" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Print a list of all flag/mode/parameter/word settings.
  

<br/><br/><strong>Main:showChangedParticleData</strong>  <input type="radio" name="13" value="on"><strong>On</strong>
<input type="radio" name="13" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Print a list of particle and decay data for those particles 
that were changed (one way or another).
  

<br/><br/><strong>Main:showAllParticleData</strong>  <input type="radio" name="14" value="on"><strong>On</strong>
<input type="radio" name="14" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Print a list of all particle and decay data.
  

<br/><br/><strong>Main:writeChangedSettings</strong>  <input type="radio" name="15" value="on"><strong>On</strong>
<input type="radio" name="15" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Write a file with the changed flag/mode/parameter/word settings, in
a format appropriate to be read in at the beginning of a new  
run, using the <code>pythia.readFile("fileName")</code> method. 
  

<br/><br/><table><tr><td><strong>Main:changedSettingsFile  </td><td></td><td> <input type="text" name="16" value="currentSettings.cmnd" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>currentSettings.cmnd</strong></code>)</td></tr></table>
The name of the file to which the changed flag/mode/parameter/word
settings are written if <code>Main:writeChangedSettings</code>
is on. 
  

<br/><br/><strong>Main:writeAllSettings</strong>  <input type="radio" name="17" value="on"><strong>On</strong>
<input type="radio" name="17" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Write a file with all flag/mode/parameter/word settings, in
a format appropriate to be read in at the beginning of a new  
run, using the <code>pythia.readFile("fileName")</code> method. 
  

<br/><br/><table><tr><td><strong>Main:allSettingsFile  </td><td></td><td> <input type="text" name="18" value="allSettings.cmnd" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>allSettings.cmnd</strong></code>)</td></tr></table>
The name of the file to which a flag/mode/parameter/word 
settings are written if <code>Main:writeAllSettings</code>
is on. 
  

<br/><br/><strong>Main:showAllStatistics</strong>  <input type="radio" name="19" value="on"><strong>On</strong>
<input type="radio" name="19" value="off" checked="checked"><strong>Off</strong>
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

<li><code>main06.cc</code> : tests of internally implemented 
cross sections for elastic and diffractive topologies, using 
<code>main06.cmnd</code> to pick process.</li>

<li><code>main07.cc</code> : tests of internally implemented 
cross sections for minimum-bias events, using 
<code>main07.cmnd</code> to pick options.</li>

<li><code>main08.cc</code> : generation of the QCD jet cross section
by splitting the run into subruns, each in its own <i>pT</i> bin,
and adding the results properly reweighted.</li>

<li><code>main09.cc</code> : generation of LEP1 hadronic events, i.e. 
<i>e^+e^- -> gamma*/Z^0 -> q qbar</i>, with charged multiplicity, 
sphericity, thrust and jet analysis. </li>

<li><code>main10.cc</code> : illustration how userHooks can be used
interact directly with the event-generation process.</li>

<li><code>main11.cc</code> : a simple example how the Les Houches
Accord interface, plus a few more Fortran-to-C++ commands, allows
hard processes to be generated by PYTHIA 6.4 and then processed 
further by PYTHIA 8.</li>

<li><code>main12.cc</code> : a fairly extensive study of 
event properties, with hard processes generated by PYTHIA 6.4. 
It reads in a <code>main12.fcmnd</code> file with commands specfically
for the Fortran PYTHIA 6.4 program and another <code>main12.ccmnd</code> 
file illustrating several of the settings listed on these pages.</li>

<li><code>main13.f</code> : a Fortran program (!) showing how 
PYTHIA 6.4 can be used to generate a Les Houches Event File 
<code>ttbar.lhe</code> with top events. This program can easily be 
modified to generate other files, bigger and/or for other processes.</li>

<li><code>main14.cc</code> : a study of top events, fed in from the 
Les Houches Event File generated by <code>main13.f</code>. This file 
currently only contains 100 events so as not to make the distributed
PYTHIA package too big, and so serves mainly as a demonstration of 
the principles involved. </li> 

<li><code>main15.cc</code> : an example how the Les Houches Accord 
interface can be used to input various toy parton-level configurations,
e.g. to study the hadronization of junction topologies.</li>

<li><code>main16.cc</code> : tests of internally implemented cross sections
for Supersymmetric particle production, with SYSY spectrum defined in
<code>main16.spc</code> and settings in <code>main16.cmnd</code>.</li>

<li><code>main17.cc</code> : shows how an external decay handler can 
be linked to handle the decays of some particles.</li>

<li><code>main18.cc</code> : shows how an external random number 
generator can be linked to handle this task.</li>

<li><code>main21.cc</code> : similar to main01, except that the 
event record is output in the HepMC event record format. Requires 
that HepMC and CLHEP are properly linked.</li>

<li><code>main22.cc</code> : similar to main12, except that the 
event record is output in the HepMC event record format. Requires 
that PYTHIA 6.4, HepMC and CLHEP are properly linked.</li>

<li><code>main23.cc</code> : a streamlined version, where all input on
the event generation is taken from the <code>main23.cmnd</code> file,
and the ony thing done is to generate events and store them in a
file in the HepMC event record format. All physics studies will have
to be done afterwards. Requires that HepMC and CLHEP are properly 
linked.</li>

<li><code>main31.cc</code> : a test of the shape of parton densities,
as a check prior to using a given PDF set in a generator.  Requires 
that LHAPDF is properly linked.</li>

<li><code>main32.cc</code> : compares the charged multiplicity 
distribution, and a few other aspects, between default PYTHIA PDF and 
another one. Requires that LHAPDF is properly linked.</li>

</ul>

<h3>Spares</h3>

For currently unforeseen purposes, a few dummy settings are made 
available here. The user can set the desired value in a "cards file"
and then use that value in the main program as desired.

<br/><br/><strong>Main:spareFlag1</strong>  <input type="radio" name="20" value="on"><strong>On</strong>
<input type="radio" name="20" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
  

<br/><br/><strong>Main:spareFlag2</strong>  <input type="radio" name="21" value="on"><strong>On</strong>
<input type="radio" name="21" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
  

<br/><br/><strong>Main:spareFlag3</strong>  <input type="radio" name="22" value="on"><strong>On</strong>
<input type="radio" name="22" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
  

<br/><br/><table><tr><td><strong>Main:spareMode1  </td><td></td><td> <input type="text" name="23" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>)</td></tr></table>
</modeopen>

<br/><br/><table><tr><td><strong>Main:spareMode2  </td><td></td><td> <input type="text" name="24" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>)</td></tr></table>
</modeopen>

<br/><br/><table><tr><td><strong>Main:spareMode3  </td><td></td><td> <input type="text" name="25" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>)</td></tr></table>
</modeopen>

<br/><br/><table><tr><td><strong>Main:spareParm1 </td><td></td><td> <input type="text" name="26" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
  

<br/><br/><table><tr><td><strong>Main:spareParm2 </td><td></td><td> <input type="text" name="27" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
  

<br/><br/><table><tr><td><strong>Main:spareParm3 </td><td></td><td> <input type="text" name="28" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
  

<br/><br/><table><tr><td><strong>Main:spareWord1  </td><td></td><td> <input type="text" name="29" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
  

<br/><br/><table><tr><td><strong>Main:spareWord2  </td><td></td><td> <input type="text" name="30" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
  

<br/><br/><table><tr><td><strong>Main:spareWord3  </td><td></td><td> <input type="text" name="31" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
  

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
if($_POST["4"] != "1960.")
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
if($_POST["7"] != "1000")
{
$data = "Main:numberOfEvents = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "2")
{
$data = "Main:numberToList = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "50")
{
$data = "Main:timesToShow = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "10")
{
$data = "Main:timesAllowErrors = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "on")
{
$data = "Main:showChangedSettings = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "off")
{
$data = "Main:showAllSettings = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "off")
{
$data = "Main:showChangedParticleData = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "off")
{
$data = "Main:showAllParticleData = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "off")
{
$data = "Main:writeChangedSettings = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "currentSettings.cmnd")
{
$data = "Main:changedSettingsFile = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "off")
{
$data = "Main:writeAllSettings = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "allSettings.cmnd")
{
$data = "Main:allSettingsFile = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "off")
{
$data = "Main:showAllStatistics = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "off")
{
$data = "Main:spareFlag1 = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "off")
{
$data = "Main:spareFlag2 = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "off")
{
$data = "Main:spareFlag3 = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "0")
{
$data = "Main:spareMode1 = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "0")
{
$data = "Main:spareMode2 = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "0")
{
$data = "Main:spareMode3 = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
if($_POST["26"] != "0.")
{
$data = "Main:spareParm1 = ".$_POST["26"]."\n";
fwrite($handle,$data);
}
if($_POST["27"] != "0.")
{
$data = "Main:spareParm2 = ".$_POST["27"]."\n";
fwrite($handle,$data);
}
if($_POST["28"] != "0.")
{
$data = "Main:spareParm3 = ".$_POST["28"]."\n";
fwrite($handle,$data);
}
if($_POST["29"] != "void")
{
$data = "Main:spareWord1 = ".$_POST["29"]."\n";
fwrite($handle,$data);
}
if($_POST["30"] != "void")
{
$data = "Main:spareWord2 = ".$_POST["30"]."\n";
fwrite($handle,$data);
}
if($_POST["31"] != "void")
{
$data = "Main:spareWord3 = ".$_POST["31"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright C 2007 Torbjorn Sjostrand -->
