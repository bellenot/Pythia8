<html>
<head>
<title>Generic Settings</title>
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

<form method='post' action='GenericSettings.php'>

<h2>Generic Settings</h2>

Here are collected a few scattered minor topics, controlled from the 
main <code>Pythia</code> class. Currently this concerns how to set the 
seed of the random numbers and how to regulate error checks on the 
generated event. 

<h3>Random numbers</h3>

The seed of the random number generator can be set here:

<br/><br/><strong>Pythia:setSeed</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Indicates whether a user-set seed should be used every time the 
<code>Pythia::init</code> routine is called. If off, the random number 
generator is initialized with its default seed at the beginning 
of the run, and never again. If on, each new <code>Pythia::init</code> 
call (should several be made in the same run) results in the random 
number being re-initialized, thereby possibly starting over with the 
same sequence, if you do not watch out.
  

<br/><br/><table><tr><td><strong>Pythia:seed  </td><td></td><td> <input type="text" name="2" value="-1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1</strong></code>; <code>maximum = 900000000</code>)</td></tr></table>
The seed to be used, if <code>setSeed</code> is on.<br/>
A negative value gives the default seed,<br/> 
a value 0 gives a random seed based on the time, and<br/>
a value between 1 and 900,000,000 a unique different random number 
sequence.
</modeopen>

<h3>Error Checks</h3>

There is also a few settings related to error checking.

<br/><br/><strong>Pythia:checkEvent</strong>  <input type="radio" name="3" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="3" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
When an event has been successfully generated, check that the 
final event record does not contain any unphysical particles, or 
nonconserved charge or energy-momentum. If this check fails, then
<code>pythia.next()</code> obtains the value <code>false</code>.
  

<br/><br/><table><tr><td><strong>Pythia:nErrList  </td><td></td><td> <input type="text" name="4" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>)</td></tr></table>
The number of erroneous events, in the above check, for which 
event listing and other detailed information will be printed. 
After that, only the normal error messages will be issued. 
Error counters are always updated, and accumulated numbers can be   
shown with <code>pythia.statistics()</code> at the end of the run.
</modeopen>

<br/><br/><table><tr><td><strong>Pythia:epTolerance </td><td></td><td> <input type="text" name="5" value="1e-5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1e-5</strong></code>)</td></tr></table>
Maximum allowed summed deviation of <i>E</i>, <i>p_x</i>, 
<i>p_y</i> and <i>p_z</i> between the incoming beams and the 
final state, as a fraction of the initial energy. (Unfortunetely 
roundoff errors do not scale linearly with the energy, and also have 
a very long tail. So while most events at lower energies may be correct 
to better than 1e-10, at LHC it does not have to signal any fundamental 
bug if also the default tolerance above is violated occasionally.)
  

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

if($_POST["1"] != "off")
{
$data = "Pythia:setSeed = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "-1")
{
$data = "Pythia:seed = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "on")
{
$data = "Pythia:checkEvent = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "0")
{
$data = "Pythia:nErrList = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "1e-5")
{
$data = "Pythia:epTolerance = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright C 2007 Torbjorn Sjostrand -->
