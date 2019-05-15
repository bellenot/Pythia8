<html>
<head>
<title>SUSY Processes</title>
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

<form method='post' action='SUSYProcesses.php'>

<h2>SUSY Processes</h2>

Here is collected all processes involving supersymmetric particle 
production, with the exception of the (extended) Higgs sector.
Since the number of separate but closely related processes is so big,
there are not switches for each separate process but only for a
reasonable set of subgroups.

<p/>
The processes in this set are not yet fully handled.

<br/><br/><b>Important note:</b> 
In order to simulate SUSY processes it is required to read in the 
couplings and masses relevant for the scenario to be studied. This 
is done with the help of the SUSY Les Houches Accord (SLHA). The 
reading of a relevant SLHA file <b>must</b> be set up, as described 
on <?php $filepath = $_GET["filepath"];
echo "<a href='SUSYLesHouchesAccord.php?filepath=".$filepath."' target='page'>";?>this page</a>.

<br/><br/><strong>SUSY:all</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for production of supersymmetric particles, i.e. 
particles with R-parity -1. 
  

<br/><br/><strong>SUSY:qqbar2chi0chi0</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production of neutralinos by quark-antiquark annihilation. With
four neutralino species this gives ten separate processes, codes 1201
- 1210. The current implementation is in a development stage and does
not include contributions from diagrams with intermediate squarks,
hence the correct answer is only obtained in the limit of infinitely heavy
squarks. As yet, there are also no sparticle decays implemented - 
all neutralinos decay to Gravitino+photon.
  


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
$data = "SUSY:all = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "SUSY:qqbar2chi0chi0 = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2007 Torbjorn Sjostrand -->

