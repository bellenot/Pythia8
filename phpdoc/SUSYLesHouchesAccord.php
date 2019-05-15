<html>
<head>
<title>SUSY Les Houches Accord</title>
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

<form method='post' action='SUSYLesHouchesAccord.php'>

<h2>SLHA</h2>

The PYTHIA 8 program does not contain an internal spectrum calculator
(a.k.a. RGE package) to provide supersymmetric couplings, mixing angles,
masses and branching ratios. Thus the SUSY Les Houches Accord (SLHA)
[<a href="Bibliography.php" target="page">Ska04</a>][<a href="Bibliography.php" target="page">All08</a>] is the only way of
inputting SUSY models, and SUSY processes cannot be run unless such an
input has taken place. 

Most of the SUSY implementation in PYTHIA 8 is compatible with both the 
SLHA1 and SLHA2 conventions 
(with the exception of R-parity violation and the NMSSM
extension in the latter case). Internally, PYTHIA 8 uses the 
SLHA2 conventions and translates SLHA1 input to these when necessary. 
See the section on SUSY Processes for more information.

When reading LHEF files, Pythia automatically looks for SLHA information
between <code>&lt;slha&gt;...&lt;/slha&gt;</code> tags in the header of such
files. When running Pythia without LHEF input (or if reading an LHEF
file that does not contain SLHA information in the header), 
a separate file containing SLHA information may be specified
using <code>SLHA:file</code> (see below).  

With the so-called <code>QNUMBERS</code>
extension [<a href="Bibliography.php" target="page">Alw07</a>], the SLHA input format can also be used for
more general BSM models, although the implementation of this extension is not
yet complete in PYTHIA 8. 

Finally, the SLHA input capability can of
course also be used to input SLHA-formatted MASS and DECAY tables for other
particles, such as the top quark, furnishing a less sophisticated but more
universal complement to the
standard PYTHIA 8-specific methods for inputting such information (for the
latter, see the section on Particle Data). 

<h3>SLHA Switches and Parameters</h3>

<br/><br/><table><tr><td><strong>SLHA:file  </td><td></td><td> <input type="text" name="1" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
Name of an SLHA (or LHEF) file containing the
SUSY/BSM model definition, spectra, and (optionally) decay tables.
  

<br/><br/><table><tr><td><strong>SLHA:minMassSM </td><td></td><td> <input type="text" name="2" value="100.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>100.0</strong></code>)</td></tr></table>
Since the default output of some programs includes SLHA output for
particles that are not desired changed in Pythia, such as the c quark,
this parameter provides a possibility to ignore SLHA input for all SM
particles (speficially particles with PDG codes below 1 million) 
whose default masses in Pythia lie below some threshold value, given
by this parameter. The default value of 100.0 allows SLHA input to
modify the top quark, but not, e.g., the Z and W bosons. 
  

<p/><code>mode&nbsp; </code><strong> SLHA:verbose &nbsp;</strong> 
 (<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 3</code>)<br/>
Controls amount of text output written by the SLHA interface, with a
value of 0 corresponding to the most quiet mode.
  

<h2>Internal SLHA Variables</h2>
The following variables are used internally by PYTHIA as local copies
of SLHA information. User changes will generally have no effect, since
these variables will be reset by the SLHA reader during initialization.
<br/><br/><strong>SLHA:NMSSM</strong>  <input type="radio" name="3" value="on"><strong>On</strong>
<input type="radio" name="3" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Corresponds to SLHA block MODSEL entry 3.
  


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

if($_POST["1"] != "void")
{
$data = "SLHA:file = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "100.0")
{
$data = "SLHA:minMassSM = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "off")
{
$data = "SLHA:NMSSM = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2009 Torbjorn Sjostrand -->


