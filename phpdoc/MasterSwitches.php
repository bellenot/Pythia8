<html>
<head>
<title>Master Switches</title>
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

<form method='post' action='MasterSwitches.php'>

<h2>Master Switches</h2>

Sometimes it may be convenient to omit certain aspects of the event 
generation chain. This cannot be motivated in a full-blown production
run, but can often be convenient for own understanding and for
debug purposes. The flags on this page allow just that.

<p/>
The event generation is subdivided into three levels: the process
level, the parton level and the hadron level, and flags are grouped
accordingly. 

<h3>Process Level</h3>

The <code>ProcessLevel</code> class administrates the initial step of 
the event generation, wherein the basic process is selected. Currently 
this is done either using some of the internal processes, or with 
Les Houches Accord input.

Since there cannot be any event at all without an initial process,
there is no possibility to switch off this part of the story. It is
possible, however, to stop the generation immediately after the
basic process has been selected:

<br/><br/><strong>Pythia:partonLevel</strong>  <input type="radio" name="1" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="1" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
If off then stop the generation after the hard process has been 
generated, but before the parton-level and hadron-level steps. 
The <code>process</code> record is filled, but the <code>event</code> 
one is not.
  
 
<h3>PartonLevel</h3>

The <code>PartonLevel</code> class administrates the middle step of the 
event generation, i.e. the evolution from an input (hard) process from
<code>ProcessLevel</code>, containing a few partons only, to a complete 
parton-level configuration to be handed on to <code>HadronLevel</code>. 
This step involves the application of initial- and final-state radiation, 
multiple interactions and the structure of beam remnants.

<p/>
Some parts of the event generation on this level may be switched off
individually: 

<br/><br/><strong>PartonLevel:MI</strong>  <input type="radio" name="2" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="2" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Master switch for multiple interactions; on/off = true/false.
Further options are found <?php $filepath = $_GET["filepath"];
echo "<a href='MultipleInteractions.php?filepath=".$filepath."' target='page'>";?>here</a>.
  

<br/><br/><strong>PartonLevel:ISR</strong>  <input type="radio" name="3" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="3" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Master switch for initial-state radiation; on/off = true/false.
Further options are found <?php $filepath = $_GET["filepath"];
echo "<a href='SpacelikeShowers.php?filepath=".$filepath."' target='page'>";?>here</a>.
  

<br/><br/><strong>PartonLevel:FSR</strong>  <input type="radio" name="4" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="4" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Master switch for initial-state radiation; on/off = true/false.
Further options are found <?php $filepath = $_GET["filepath"];
echo "<a href='TimelikeShowers.php?filepath=".$filepath."' target='page'>";?>here</a>.
If you leave this switch on, the following two switches allow 
more detailed control to switch off only parts of the showers. 
  

<br/><br/><strong>PartonLevel:FSRinProcess</strong>  <input type="radio" name="5" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="5" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Switch for final-state radiation in association with the hard process 
itself; on/off = true/false. In addition <code>PartonLevel:FSR</code>
must be on for these emissions to occur. 
  

<br/><br/><strong>PartonLevel:FSRinResonances</strong>  <input type="radio" name="6" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="6" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Master switch for final-state radiation in any resonance decays 
subsequent to the hard process itself; on/off = true/false. In addition 
<code>PartonLevel:FSR</code> must be on for these emissions to occur.
  

<p/>
Switching off all the above switches is <i>not</i> equivalent to 
setting <code>Pythia:partonLevel = off</code>. In the former case a minimal 
skeleton of parton-level operations are carried out, such as tying together 
the scattered partons with the beam remnants into colour singlets, and
storing this information in the <code>event</code> record. It is therefore
possible to go on and hadronize the event, if desired. In the latter case
<i>no</i> operations at all are carried out on the parton level, and 
therefore it is also not possible to go on to the hadron level.

<p/>
Finally, it is possible to stop the generation immediately after the
parton-level step:

<br/><br/><strong>Pythia:hadronLevel</strong>  <input type="radio" name="7" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="7" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
If off then stop the generation after the hard process and 
parton-level activity has been generated, but before the 
hadron-level steps.
  

<h3>HadronLevel</h3>

The <code>HadronLevel</code> class administrates the final step of the 
event generation, wherein the partonic configuration from 
<code>PartonLevel</code> is hadronized, including string fragmentation 
and secondary decays.

<p/>
Most of the code in this class deals with subdividing the partonic
content of the event into separate colour singlets, that can be
treated individually by the string fragmentation machinery. When a
junction and an antijunction are directly connected, it also breaks 
the string between the two, so that the topology can be reduced back 
to two separate one-junction systems, while still preserving the
expected particle flow in the junction-junction string region(s).

<p/>
Some parts of the event generation on this level may be switched off
individually: 

<br/><br/><strong>HadronLevel:Hadronize</strong>  <input type="radio" name="8" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="8" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Master switch for hadronization; on/off = true/false.
Further options are found <?php $filepath = $_GET["filepath"];
echo "<a href='Fragmentation.php?filepath=".$filepath."' target='page'>";?>here</a>.
  

<br/><br/><strong>HadronLevel:Decay</strong>  <input type="radio" name="9" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="9" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Master switch for decays; on/off = true/false.
Further options are found <?php $filepath = $_GET["filepath"];
echo "<a href='ParticleDecays.php?filepath=".$filepath."' target='page'>";?>here</a>.
  

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

if($_POST["1"] != "on")
{
$data = "Pythia:partonLevel = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "on")
{
$data = "PartonLevel:MI = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "on")
{
$data = "PartonLevel:ISR = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "on")
{
$data = "PartonLevel:FSR = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "on")
{
$data = "PartonLevel:FSRinProcess = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "on")
{
$data = "PartonLevel:FSRinResonances = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "on")
{
$data = "Pythia:hadronLevel = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "on")
{
$data = "HadronLevel:Hadronize = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "on")
{
$data = "HadronLevel:Decay = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright C 2007 Torbjorn Sjostrand -->
