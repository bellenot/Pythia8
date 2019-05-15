<html>
<head>
<title>SUSY Processes</title>
<link rel="stylesheet" type="text/css" href="pythia.css"/>
<link rel="shortcut icon" href="pythia32.gif"/>
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

The implementation of SUSY processes is barely begun, so this page
is more a placeholder than a repository of useful processes.

<p/>
Here is collected processes involving supersymmetric particle 
production, with the exception of the (extended) Higgs sector.
Since the number of separate but closely related processes is so big,
there will not be switches for each separate process but only for a
reasonable set of subgroups.

<br/><br/><b>Important note:</b> 
In order to simulate SUSY processes it is required to read in the 
couplings and masses relevant for the scenario to be studied. This 
is done with the help of the SUSY Les Houches Accord (SLHA), including
the SLHA2 extensions and generalizations. (Internally, the SLHA2
conventions are used. SLHA1 spectra are automatically translated into
SLHA2 notation during initialization.) The 
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
four neutralino species this gives ten separate processes, codes 
1201 - 1210. The cross section expressions follow [<a href="Bibliography.php" target="page">Boz07</a>] and
should thus be  valid also in the case of non-minimal flavour
violation and/or CP violation. 
  

<br/><br/><strong>SUSY:qqbar2chi+-chi0</strong>  <input type="radio" name="3" value="on"><strong>On</strong>
<input type="radio" name="3" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Associated chargino-neutralino production by quark-antiquark
annihilation. With four neutralino species, two chargino ones, and
maintaining charge conjugate proceeses separate, this gives 16 
separate processes, codes 1221 - 1236. The cross section expressions 
follow [<a href="Bibliography.php" target="page">Boz07</a>] and should thus be valid also in the case of 
non-minimal flavour violation and/or CP violation. 
  

<br/><br/><strong>SUSY:qqbar2chi+chi-</strong>  <input type="radio" name="4" value="on"><strong>On</strong>
<input type="radio" name="4" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production of charginos by quark-antiquark annihilation. With
two chargino species and maintaining mutually charge conjugate
processes separate, this gives four separate processes, codes 
1241 - 1244. The cross section expressions follow [<a href="Bibliography.php" target="page">Boz07</a>] 
and should thus be valid also in the case of non-minimal flavour 
violation and/or CP violation. 
  

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
if($_POST["3"] != "off")
{
$data = "SUSY:qqbar2chi+-chi0 = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "off")
{
$data = "SUSY:qqbar2chi+chi- = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2009 Peter Skands, Torbjorn Sjostrand -->

