<html>
<head>
<title>Higgs Processes</title>
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

<form method='post' action='HiggsProcesses.php'>

<h2>Higgs Processes</h2>

This page contains Higgs production within and beyond the Standard Model.
It is at a very early stage.

<h3>Standard-Model Higgs</h3>

<br/><br/><strong>SMHiggs:all</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of Higgs production within the Standard Model.
  

<br/><br/><strong>SMHiggs:ffbar2H</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> H^0</i>, where <i>f</i> sums over available
flavours except top.
Code 801.
  

<br/><br/><strong>SMHiggs:gg2H</strong>  <input type="radio" name="3" value="on"><strong>On</strong>
<input type="radio" name="3" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> H^0</i>.
Code 802.
  

<br/><br/><strong>SMHiggs:gmgm2H</strong>  <input type="radio" name="4" value="on"><strong>On</strong>
<input type="radio" name="4" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>gamma gamma -> H^0</i>.
Code 803.
  

<br/><br/><strong>SMHiggs:ffbar2HZ</strong>  <input type="radio" name="5" value="on"><strong>On</strong>
<input type="radio" name="5" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> H^0 Z^0</i>.
Code 804.
  

<br/><br/><strong>SMHiggs:ffbar2HW</strong>  <input type="radio" name="6" value="on"><strong>On</strong>
<input type="radio" name="6" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> H^0 W^+-</i>.
Code 805.
  

<br/><br/><strong>SMHiggs:qg2Hq</strong>  <input type="radio" name="7" value="on"><strong>On</strong>
<input type="radio" name="7" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q g -> H^0 q</i>. This process gives first-order 
corrections to the <i>f fbar -> H^0</i> one above (for 
<i>f = q </i>), so both cannot be used simultaneously without 
unphysical doublecounting. The current one should only be used to study 
the high-pT tail, while <i>f g -> H^0 f</i> should be used for 
inclusive production. Therefore this process is not switched on by 
the <code>SMHiggs:all</code> flag.
Code 806.
  

<h3>Parameters for Higgs production and decay</h3>

<br/><br/><strong>ResonanceH::linearWidthWWZZ</strong>  <input type="radio" name="8" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="8" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
The partial width of a Higgs particle to a pair of gauge bosons,
<i>W^+ W^-</i> or <i>Z^0 Z^0</i>, depends cubically on the
Higgs mass. When selecting the Higgs according to a Breit-Wigner,
so that the actual mass <i>mHat</i> does not agree with the
nominal <i>m_Higgs</i> one, an ambiguity arises which of the 
two to use [<a href="Bibliography.php" target="page">Sey95</a>]. The default is to use a linear 
dependence on <i>mHat</i>, i.e. a width proportional to 
<i>m_Higgs^2 * mHat</i>, while <code>off</code> gives a 
<i>mHat^3</i> dependence. This does not affect the width to 
fermions, which only depends on <i>mHat</i>.
  

<br/><br/><table><tr><td><strong>SMHiggs:parity  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 3</code>)</td></tr></table>
<modepick name="SMHiggs:parity" default="1" min="0" max="3">
possibility to modify angular decay correlations in Higgs decay to
<ei>Z^0 Z^0</ei> or <ei>W^+ W^-</ei> to four fermions. Currently it 
does not affect the partial width of the channels.
<br/>
<input type="radio" name="9" value="0"><strong>0 </strong>: isotropic decays.<br/>
<input type="radio" name="9" value="1" checked="checked"><strong>1 </strong>: assuming the Higgs is a pure scalar (CP-even), as in the Standard Model.<br/>
<input type="radio" name="9" value="2"><strong>2 </strong>: assuming the Higgs is a pure pseudoscalar(CP-odd).<br/>
<input type="radio" name="9" value="3"><strong>3 </strong>: assuming the Higgs is a mixture of the two, including the CP-violating interference term. The parameter<ei>eta</ei>, see below, sets the strength of the CP-odd admixture,with the interference term being proportional to <ei>eta</ei>and the CP-odd one to <ei>eta^2</ei>.<br/>
</modepick>

<br/><br/><table><tr><td><strong>SMHiggs:etaParity </td><td></td><td> <input type="text" name="10" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
The <i>eta</i> value of CP-violation in the 
<code>SMHiggs:parity = 2</code> option. 
  

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
$data = "SMHiggs:all = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "SMHiggs:ffbar2H = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "off")
{
$data = "SMHiggs:gg2H = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "off")
{
$data = "SMHiggs:gmgm2H = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "off")
{
$data = "SMHiggs:ffbar2HZ = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "off")
{
$data = "SMHiggs:ffbar2HW = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "off")
{
$data = "SMHiggs:qg2Hq = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "on")
{
$data = "ResonanceH::linearWidthWWZZ = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "1")
{
$data = "SMHiggs:parity = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "0.")
{
$data = "SMHiggs:etaParity = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright C 2007 Torbjorn Sjostrand -->

