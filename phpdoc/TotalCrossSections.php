<html>
<head>
<title>Total Cross Sections</title>
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

<form method='post' action='TotalCrossSections.php'>

<h2>Total Cross Sections</h2>

The <code>SigmaTotal</code> class returns the total, elastic, diffractive 
and nondiffractive cross sections in hadronic collisions, and also the
slopes of the <i>d(sigma)/dt</i> distributions. The parametrizations 
used are from [<a href="Bibliography.php" target="page">Sch97</a>] which borrows some of the total cross 
sections from [<a href="Bibliography.php" target="page">Don92</a>].

<p/>
The allowed combinations of incoming particles are <i>p + p</i>, 
<i>pbar + p</i>, <i>pi+ + p</i>, <i>pi- + p</i>, 
<i>pi0/rho0 + p</i>, <i>phi + p</i>, <i>J/psi + p</i>, 
<i>rho + rho</i>, <i>rho + phi</i>, <i>rho + J/psi</i>, 
<i>phi + phi</i>, <i>phi + J/psi</i>, <i>J/psi + J/psi</i>.   
The strong emphasis on vector mesons is related to the description
of <i>gamma + p</i> and <i>gamma + gamma</i> interactions in a 
Vector Dominance Model framework (which will not be available for some 
time to come, so this is a bit of overkill).

<h3>Variables</h3>

If the internally implemented cross section parametrizations are not 
satisfactory, it is possible to override the cross section values 
(but currently not the <i>t</i> slopes), with 

<br/><br/><strong>SigmaTotal:setOwn</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>no</strong></code>)<br/>
Allow a user to set own cross sections by hand; yes/no = true/false.
  

<p/>
When <code>SigmaTotal:setOwn = yes</code>, the user is expected to set 
values for the corresponding cross sections:

<br/><br/><table><tr><td><strong>SigmaTotal:sigmaTot </td><td></td><td> <input type="text" name="2" value="80." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>80.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
Total cross section in mb.
  

<br/><br/><table><tr><td><strong>SigmaTotal:sigmaEl </td><td></td><td> <input type="text" name="3" value="20." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>20.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
Elastic cross section in mb.
  

<br/><br/><table><tr><td><strong>SigmaTotal:sigmaXB </td><td></td><td> <input type="text" name="4" value="8." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>8.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
Single Diffractive cross section <i>A + B -> X + B</i> in mb.
  

<br/><br/><table><tr><td><strong>SigmaTotal:sigmaAX </td><td></td><td> <input type="text" name="5" value="8." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>8.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
Single Diffractive cross section <i>A + B -> A + X</i> in mb.
  

<br/><br/><table><tr><td><strong>SigmaTotal:sigmaXX </td><td></td><td> <input type="text" name="6" value="4." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>4.</strong></code>; <code>minimum = 0.</code>)</td></tr></table>
Double Diffractive cross section <i>A + B -> X_1 + X_2</i> in mb.
  

<p/>
Note that the total cross section subtracted by the elastic and various 
diffractive ones gives the inelastic nondiffractive cross section, 
which therefore is not set separately. If this cross section evaluates 
to be negative the internal parametrizations are used instead of the 
ones here. However, since the nondiffractive inelastic cross section 
is what makes up the minimum-bias event class, and plays a major role 
in the description of multiple interactions, it is important that a 
consistent set is used. 

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

if($_POST["1"] != "no")
{
$data = "SigmaTotal:setOwn = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "80.")
{
$data = "SigmaTotal:sigmaTot = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "20.")
{
$data = "SigmaTotal:sigmaEl = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "8.")
{
$data = "SigmaTotal:sigmaXB = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "8.")
{
$data = "SigmaTotal:sigmaAX = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "4.")
{
$data = "SigmaTotal:sigmaXX = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2007 Torbjorn Sjostrand -->
