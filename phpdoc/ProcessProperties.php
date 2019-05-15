<html>
<head>
<title>Process Properties</title>
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

<form method='post' action='ProcessProperties.php'>

<h2>Process Properties</h2>

Here is collected some possibilities to affect the generation of 
all internally implemented processes in one go, or some coherent
subset thereof. Phase space cuts appear on a special page.

<h3>Incoming partons</h3>

<code>InFlux</code> is base class for the combination of allowed incoming 
partons in a given process, and keeps track of parton densities and 
common weight factors (such as charge) required to give the process 
cross section.

<p/>
There is one useful degree of freedom in this class:
<br/><br/><table><tr><td><strong>InFlux:nQuark  </td><td></td><td> <input type="text" name="1" value="5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>5</strong></code>; <code>minimum = 0</code>; <code>maximum = 5</code>)</td></tr></table>
Number of allowed incoming quark flavours in the beams; a change 
to 4 would thus exclude <i>b</i> and <i>bbar</i> as incoming 
partons, etc.
</modeopen>

<p/>
There is also the possibility to obtain some documentation.
<br/><br/><strong>InFlux:showChannels</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
At initialization show which incoming flavours and flavour combinations
have beeen set up, process by process, along with initial weight 
(e.g. charges or CKM factors) assigned to each of them.
  

<h3>Generic cross sections</h3>

<code>SigmaProcess</code> is base class for all hard processes 
implemented in PYTHIA 8. For <i>2 -> 1</i> processes it 
should give <i>sigmaHat(sHat)</i>, for <i>2 -> 2</i> ones 
<i>d(sigmaHat(sHat, tHat))/d(tHat)</i>. The matrix-element
coding is also used by the multiple-interactions machinery,
but with a separate choice of <i>alpha_strong(M_Z^2)</i> value
and running, and separate scale choices. Also, in <i>2 -> 2</i> 
processes where resonances are produced, their couplings and thereby 
their Breit-Wigner shapes are always evaluated with the resonance mass 
as scale, irrespective of the choices below.

<p/>
The size of QCD cross sections is mainly determined by 
<br/><br/><table><tr><td><strong>SigmaProcess:alphaSvalue </td><td></td><td> <input type="text" name="3" value="0.1265" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.1265</strong></code>; <code>minimum = 0.06</code>; <code>maximum = 0.25</code>)</td></tr></table>
The <i>alpha_strong</i> value at scale <i>M_Z^2</i>. 
  

<p/>
The actual value is then regulated by the running to the <i>Q^2</i> 
renormalization scale, at which <i>alpha_strong</i> is evaluated
<br/><br/><table><tr><td><strong>SigmaProcess:alphaSorder  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
<modepick name="SigmaProcess:alphaSorder" default="1" min="0" max="2">
Order at which <ei>alpha_strong</ei> runs,
<br/>
<input type="radio" name="4" value="0"><strong>0 </strong>: zeroth order, i.e. <ei>alpha_strong</ei> is kept fixed.<br/>
<input type="radio" name="4" value="1" checked="checked"><strong>1 </strong>: first order, which is the normal value.<br/>
<input type="radio" name="4" value="2"><strong>2 </strong>: second order. Since other parts of the code do not go to second order there is no strong reason to use this option, but there is also nothing wrong with it.<br/>
</modepick>

<p/>
QED interactions are regulated by the <i>alpha_electromagnetic</i>
value at the <i>Q^2</i> renormalization scale of an interaction. 
<br/><br/><table><tr><td><strong>SigmaProcess:alphaEMorder  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = -1</code>; <code>maximum = 1</code>)</td></tr></table>
<modepick name="SigmaProcess:alphaEMorder" default="1" min="-1" max="1">
The running of <ei>alpha_em</ei> used in hard processes.
<br/>
<input type="radio" name="5" value="1" checked="checked"><strong>1 </strong>: first-order running, constrained to agree with<code>StandardModel:alphaEMmZ</code> at the <ei>Z^0</ei> mass.<br/>
<input type="radio" name="5" value="0"><strong>0 </strong>: zeroth order, i.e. <ei>alpha_em</ei> is kept fixed at its value at vanishing momentum transfer.<br/>
<input type="radio" name="5" value="-1"><strong>-1 </strong>: zeroth order, i.e. <ei>alpha_em</ei> is kept fixed, but at <code>StandardModel:alphaEMmZ</code>, i.e. its valueat the <ei>Z^0</ei> mass.<br/>
</modepick>

<p/>
The <i>Q^2</i> renormalization scale is <i>sHat</i> for 
<i>2 -> 1</i> processes while there is a choice for <i>2 -> 2</i>
ones:  
<br/><br/><table><tr><td><strong>SigmaProcess:renormScale  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 1</code>)</td></tr></table>
<modepick name="SigmaProcess:renormScale" default="0" min="0" max="1">
The <ei>Q^2</ei> renormalization scale for <ei>2 -> 2</ei> processes.
<br/>
<input type="radio" name="6" value="0" checked="checked"><strong>0 </strong>: the product of the transverse masses of the twooutgoing particles, i.e. <ei>sqrt((pT^2 + m3^2) * (pT2 + m4^2))</ei>.<br/>
<input type="radio" name="6" value="1"><strong>1 </strong>: the average of the squared transverse masses, i.e. <ei>pT^2 + 0.5 * (m3^2 + m4^2)</ei>. Useful for comparisons with PYTHIA 6.<br/>
</modepick>

<p/>
It is also possible to shift the choice by a global factor.
<br/><br/><table><tr><td><strong>SigmaProcess:renormMult </td><td></td><td> <input type="text" name="7" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 10.</code>)</td></tr></table>
The <i>Q^2</i> renormalization scale for <i>2 -> 1</i> and
<i>2 -> 2</i> processes is multiplied by this factor relative to the
scale described above. Should be use sparingly for <i>2 -> 1</i>
processes. 

<p/>
Corresponding options exist for the <i>Q^2</i> factorization scale
used as argument in PDF's. Again <i>sHat</i> is scale for <i>2 -> 1</i> 
processes while there is a choice for <i>2 -> 2</i> ones:
<br/><br/><table><tr><td><strong>SigmaProcess:factorScale  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 1</code>)</td></tr></table>
<modepick name="SigmaProcess:factorScale" default="0" min="0" max="1">
The <ei>Q^2</ei> factorization scale for <ei>2 -> 2</ei> processes.
<br/>
<input type="radio" name="8" value="0" checked="checked"><strong>0 </strong>: the smaller of the squared transverse masses of the twooutgoing particles, i.e. <ei>pT^2 + min(m3^2 + m4^2)</ei>.<br/>
<input type="radio" name="8" value="1"><strong>1 </strong>: the average of the squared transverse masses, i.e. <ei>pT^2 + 0.5 * (m3^2 + m4^2)</ei>. Useful for comparisons with PYTHIA 6.<br/>
</modepick>

<p/>
It is also possible to shift the choice by a global factor.
<br/><br/><table><tr><td><strong>SigmaProcess:factorMult </td><td></td><td> <input type="text" name="9" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>; <code>minimum = 0.1</code>; <code>maximum = 10.</code>)</td></tr></table>
The <i>Q^2</i> factorization scale for <i>2 -> 1</i> and
<i>2 -> 2</i> processes is multiplied by this factor relative to the
scale described above. Should be use sparingly for <i>2 -> 1</i>
processes.  

<h3>Special cross sections</h3>

Here settings that affect some special group of processes, but not all.

<br/><br/><table><tr><td><strong>SigmaProcess:nQuark  </td><td></td><td> <input type="text" name="10" value="3" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>3</strong></code>; <code>minimum = 0</code>; <code>maximum = 5</code>)</td></tr></table>
Number of allowed outgoing new quark flavours in 
<i>q qbar -> q' qbar'</i> and <i>g g-> q qbar</i> processes, 
where quarks are treated as massless in the matrix-element expressions 
(but correctly in the phase space). It is thus assumed that <i>c cbar</i> 
and <i>b bbar</i> are added separately with masses taken into account. 
A change to 4 would also include <i>c cbar</i> in the massless 
approximation, etc. 
</modeopen>

<br/><br/><table><tr><td><strong>SigmaProcess:nQuarkInLoop  </td><td></td><td> <input type="text" name="11" value="5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>5</strong></code>; <code>minimum = 3</code>; <code>maximum = 6</code>)</td></tr></table>
Number of quark flavours included in the box graphs resposible for 
<i>g g -> g gamma</i> and <i>g g-> gamma gamma</i> processes.
Owing to the complexity if the massive expressions, quarks are treated 
as massless. The default value should be applicable in the range of 
transverse momenta above the <i>b</i> mass but below the <i>t</i> one.
</modeopen>

<br/><br/><table><tr><td><strong>SigmaProcess:gmZmode  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
<modepick name="SigmaProcess:gmZmode" default="0" min="0" max="2">
Choice of full <ei>gamma^*/Z^0</ei> structure or not in relevant 
processes.
<br/>
<input type="radio" name="12" value="0" checked="checked"><strong>0 </strong>: full <ei>gamma^*/Z^0</ei> structure,with interference included.<br/>
<input type="radio" name="12" value="1"><strong>1 </strong>: only pure <ei>gamma^*</ei> contribution.<br/>
<input type="radio" name="12" value="2"><strong>2 </strong>: only pure <ei>Z^0</ei> contribution.<br/>
</modepick>

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

if($_POST["1"] != "5")
{
$data = "InFlux:nQuark = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "InFlux:showChannels = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "0.1265")
{
$data = "SigmaProcess:alphaSvalue = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "1")
{
$data = "SigmaProcess:alphaSorder = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "1")
{
$data = "SigmaProcess:alphaEMorder = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "0")
{
$data = "SigmaProcess:renormScale = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "1.")
{
$data = "SigmaProcess:renormMult = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "0")
{
$data = "SigmaProcess:factorScale = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "1.")
{
$data = "SigmaProcess:factorMult = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "3")
{
$data = "SigmaProcess:nQuark = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "5")
{
$data = "SigmaProcess:nQuarkInLoop = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "0")
{
$data = "SigmaProcess:gmZmode = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright C 2007 Torbjorn Sjostrand -->
