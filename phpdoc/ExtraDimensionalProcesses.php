<html>
<head>
<title>Extra-Dimensional Processes</title>
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

<form method='post' action='ExtraDimensionalProcesses.php'>

<h2>Extra-Dimensional Processes</h2>

Scenarios with extra dimensions allow a multitude of processes.
Currently two different categories of processes are implemented. 
The first involves the production of an excited graviton state 
<i>G^*</i> within a Randall-Sundrum (RS) scenario, the second 
phenomena from large extra dimensions (LED). Due to the close 
relation between the LED model and a so-called unparticle model, 
unparticle processes are also kept in this section.


<h3>Randall-Sundrum Resonances, production processes</h3>

The Graviton resonance state is assigned PDG code 5100039. Decays 
into fermion, gluon and photon pairs are handled with the correct 
angular distributions, while other decay channels currently are 
handled isotropically.

<p/>
There are two lowest-order processes that together normally 
should be sufficient for a simulation of <i>G^*</i> production. 

<br/><br/><strong>ExtraDimensionsG*:all</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of lowest-order <i>G^*</i> production
processes, i.e. the two ones below.
  

<br/><br/><strong>ExtraDimensionsG*:gg2G*</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>g g -> G^*</i>. 
Code 5001.
  

<br/><br/><strong>ExtraDimensionsG*:ffbar2G*</strong>  <input type="radio" name="3" value="on"><strong>On</strong>
<input type="radio" name="3" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> G^*</i>. 
Code 5002.
  

<p/>
In addition there are three first-order processes included. These 
are of less interest, but can be used for dedicated studies of the 
high-<i>pT</i> tail of <i>G^*</i> production. As usual, it would 
be double counting to include the lowest-order and first-order 
processes simultaneously. Therefore the latter ones are not included 
with the <code>ExtraDimensionsG*:all = on</code> option. In this set 
of processes all decay angles are assumed isotropic.

<br/><br/><strong>ExtraDimensionsG*:gg2G*g</strong>  <input type="radio" name="4" value="on"><strong>On</strong>
<input type="radio" name="4" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>g g -> G^* g</i>. 
Code 5003.
  

<br/><br/><strong>ExtraDimensionsG*:qg2G*q</strong>  <input type="radio" name="5" value="on"><strong>On</strong>
<input type="radio" name="5" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q g -> G^* q</i>. 
Code 5004.
  

<br/><br/><strong>ExtraDimensionsG*:qqbar2G*g</strong>  <input type="radio" name="6" value="on"><strong>On</strong>
<input type="radio" name="6" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q qbar -> G^* g</i>. 
Code 5005.
  

<h3>Randall-Sundrum Resonances, parameters</h3>

In the above scenario the main free parameter is the <i>G^*</i> mass,
which is set as usual. In addition there is one further parameter.

<br/><br/><table><tr><td><strong>ExtraDimensionsG*:kappaMG </td><td></td><td> <input type="text" name="7" value="0.054" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.054</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
dimensionless coupling, which enters quadratically in all partial
widths of the <i>G^*</i>. Is 
<i>kappa m_G* = sqrt(2) x_1 k / Mbar_Pl</i>,
where <i>x_1 = 3.83</i> is the first zero of the <i>J_1</i> Bessel 
function and <i>Mbar_Pl</i> is the modified Planck mass.
  

<h3>Large Extra Dimensions, production processes</h3>

The LED graviton, where the KK-modes normally are summed and do not 
give rise to phenomena individually, is assigned PDG code 5000039. 
The graviton emission and virtual graviton exchange processes uses 
the same implementation as the corresponding unparticle processes, 
see further details below.

<p/>
The graviton emission processes are implemented using the conventions 
described in [<a href="Bibliography.php" target="page">Ask09</a>]. As also discussed in [<a href="Bibliography.php" target="page">Ask09</a>], 
the underlying Breit-Wigner mass distribution has to be matched to the 
graviton mass spectrum in order to achieve a high MC efficiency. 

<p/>
The virtual graviton exchange processes uses the parameters defined 
in [<a href="Bibliography.php" target="page">Giu99</a>].

<p/>
The following lowest order graviton emission processes are available.

<br/><br/><strong>ExtraDimensionsLED:monojet</strong>  <input type="radio" name="8" value="on"><strong>On</strong>
<input type="radio" name="8" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of lowest-order <i>G jet</i> emission
processes, i.e. the three ones below.
  

<br/><br/><strong>ExtraDimensionsLED:gg2Gg</strong>  <input type="radio" name="9" value="on"><strong>On</strong>
<input type="radio" name="9" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>g g -> G g</i>. 
Code 5021.
  

<br/><br/><strong>ExtraDimensionsLED:qg2Gq</strong>  <input type="radio" name="10" value="on"><strong>On</strong>
<input type="radio" name="10" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q g -> G q</i>. 
Code 5022.
  

<br/><br/><strong>ExtraDimensionsLED:qqbar2Gg</strong>  <input type="radio" name="11" value="on"><strong>On</strong>
<input type="radio" name="11" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>q qbar -> G g</i>. 
Code 5023.
  

<br/><br/><strong>ExtraDimensionsLED:ffbar2GZ</strong>  <input type="radio" name="12" value="on"><strong>On</strong>
<input type="radio" name="12" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> G Z</i>. 
Code 5024.
  

<br/><br/><strong>ExtraDimensionsLED:ffbar2Ggamma</strong>  <input type="radio" name="13" value="on"><strong>On</strong>
<input type="radio" name="13" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> G gamma</i>. This process corresponds 
to the photon limit of the <i>G Z</i> process, as described in 
[<a href="Bibliography.php" target="page">Ask09</a>].
Code 5025.
  

<p/>
The following LED processes with virtual graviton exchange are 
available.

<br/><br/><strong>ExtraDimensionsLED:ffbar2gammagamma</strong>  <input type="radio" name="14" value="on"><strong>On</strong>
<input type="radio" name="14" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> (LED G*) -> gamma gamma</i>. If the 
graviton contribution is zero, the results corresponds to the 
SM contribution, i.e. equivalent to 
<code>PromptPhoton:ffbar2gammagamma</code>.
Code 5026.
  

<br/><br/><strong>ExtraDimensionsLED:gg2gammagamma</strong>  <input type="radio" name="15" value="on"><strong>On</strong>
<input type="radio" name="15" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>gg -> (LED G*) -> gamma gamma</i>. 
Code 5027.
  

<h3>Large Extra Dimensions, parameters</h3>

<br/><br/><strong>ExtraDimensionsLED:Trunc</strong>  <input type="radio" name="16" value="on"><strong>On</strong>
<input type="radio" name="16" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Option to truncate the graviton mass spectrum. Only 
concerns the graviton emission processes.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsLED:n  </td><td></td><td> <input type="text" name="17" value="2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2</strong></code>; <code>minimum = 1</code>)</td></tr></table>
Number of extra dimensions.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsLED:MD </td><td></td><td> <input type="text" name="18" value="2000." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2000.</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Fundamental scale of gravity in <i>D = 4 + n</i> dimensions.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsLED:LambdaT </td><td></td><td> <input type="text" name="19" value="2000." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2000.</strong></code>; <code>minimum = 0.0</code>)</td></tr></table>
Ultraviolet cutoff parameter for the virtual graviton exchange processes.
  

<h3>Unparticles, production processes</h3>

As mentioned above, the similar unparticle and graviton processes 
share the same implementations. The unparticle processes, however, 
only uses the dedicated unparticle parameters below. The unparticle 
is also assigned the PDG code 5000039 and is therefore called 
<i>Graviton</i> in the event record. Further details about the 
common implementations are given below.

<p/>
All unparticle processes follow the parameter conventions described in 
[<a href="Bibliography.php" target="page">Ask09</a>]. As also discussed in [<a href="Bibliography.php" target="page">Ask09</a>], the underlying 
Breit-Wigner mass distribution has to be matched to the graviton mass 
spectrum in order to achieve a high MC efficiency. 

<p/>
The following unparticle emission processes are available.

<br/><br/><strong>ExtraDimensionsUnpart:ffbar2UZ</strong>  <input type="radio" name="20" value="on"><strong>On</strong>
<input type="radio" name="20" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> U Z</i>.  
Code 5041.
  

<br/><br/><strong>ExtraDimensionsUnpart:ffbar2Ugamma</strong>  <input type="radio" name="21" value="on"><strong>On</strong>
<input type="radio" name="21" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> U gamma</i>. 
Code 5042.
  

<p/>
The following processes with virtual unparticle exchange are available.

<br/><br/><strong>ExtraDimensionsUnpart:ffbar2gammagamma</strong>  <input type="radio" name="22" value="on"><strong>On</strong>
<input type="radio" name="22" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>f fbar -> (U*) -> gamma gamma</i>. If the unparticle  
contribution is zero in the spin-2 case, the results corresponds to 
the SM contribution, i.e. equivalent to 
<code>PromptPhoton:ffbar2gammagamma</code>.
This is also the case if the model parameter settings are not applicable. 
Code 5043.
  

<br/><br/><strong>ExtraDimensionsUnpart:gg2gammagamma</strong>  <input type="radio" name="23" value="on"><strong>On</strong>
<input type="radio" name="23" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scatterings <i>gg -> (U*) -> gamma gamma</i>. 
Code 5044.
  

<h3>Unparticles, parameters</h3>

<br/><br/><strong>ExtraDimensionsUnpart:Trunc</strong>  <input type="radio" name="24" value="on"><strong>On</strong>
<input type="radio" name="24" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Option to truncate the unparticle mass spectrum.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsUnpart:spinU  </td><td></td><td> <input type="text" name="25" value="2" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Unparticle spin.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsUnpart:dU </td><td></td><td> <input type="text" name="26" value="1.4" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.4</strong></code>; <code>minimum = 1.0</code>)</td></tr></table>
Scale dimension parameter.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsUnpart:LambdaU </td><td></td><td> <input type="text" name="27" value="2000." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2000.</strong></code>; <code>minimum = 1.0</code>)</td></tr></table>
Unparticle renormalization scale.
  

<br/><br/><table><tr><td><strong>ExtraDimensionsUnpart:lambda </td><td></td><td> <input type="text" name="28" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>)</td></tr></table>
Unparticle coupling to the SM fields
  

<br/><br/><table><tr><td><strong>ExtraDimensionsUnpart:ratio </td><td></td><td> <input type="text" name="29" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 1.0</code>; <code>maximum = 1.0</code>)</td></tr></table>
Ratio, <i>lambda'/lambda</i>, between the two possible coupling constants 
of the spin-2 ME. <b>Warning:</b> A <i>ratio</i> value different from one 
give rise to an IR divergence which makes the event generation very slow, so 
this values is fixed to <i>ratio = 1</i> for the moment.
  

<h3>Common LED and unparticle implementation</h3>

Since the cross sections of the available LED and corresponding (spin-2) 
unparticle processes only differ by some constant factors they share the 
same implementation. The intention is that this should both minimize the 
amount of code needed as well as emphasize the rather small differences 
between the models from a phenomenological point of view.

The main documentation for the common LED/unparticle implementation together 
with the parameter conventions used is presented in [<a href="Bibliography.php" target="page">Ask09</a>]. Since 
this paper is focused on graviton/unparticle emission some complementary 
information, not covered in [<a href="Bibliography.php" target="page">Ask09</a>], is summarized below.

The spin-2 unparticle and graviton processes share the same matrix 
elements (ME). MEs were taken from the following papers: 
(monojets) [<a href="Bibliography.php" target="page">Giu99</a>], (mono-Z or -photon) [<a href="Bibliography.php" target="page">Che07</a>] 
(gammagamma) [<a href="Bibliography.php" target="page">Kum08</a>]. 
All unparticle processes (spin-0, spin-1 and spin-2) have an universal 
unparticle - SM coupling, <i>lambda</i>, implemented according to the 
effective operators in [<a href="Bibliography.php" target="page">Che07</a>]. The virtual graviton exchange 
processes are obtained from the spin-2 unparticle MEs by setting 
<i>dU = 2</i> and <i>lambda^2 * Chi = 4 pi</i>, which reproduces 
the graviton formulas in [<a href="Bibliography.php" target="page">Giu99</a>] and [<a href="Bibliography.php" target="page">Giu04</a>].

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
$data = "ExtraDimensionsG*:all = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "ExtraDimensionsG*:gg2G* = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "off")
{
$data = "ExtraDimensionsG*:ffbar2G* = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "off")
{
$data = "ExtraDimensionsG*:gg2G*g = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "off")
{
$data = "ExtraDimensionsG*:qg2G*q = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "off")
{
$data = "ExtraDimensionsG*:qqbar2G*g = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "0.054")
{
$data = "ExtraDimensionsG*:kappaMG = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "off")
{
$data = "ExtraDimensionsLED:monojet = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "off")
{
$data = "ExtraDimensionsLED:gg2Gg = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "off")
{
$data = "ExtraDimensionsLED:qg2Gq = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "off")
{
$data = "ExtraDimensionsLED:qqbar2Gg = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "off")
{
$data = "ExtraDimensionsLED:ffbar2GZ = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "off")
{
$data = "ExtraDimensionsLED:ffbar2Ggamma = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "off")
{
$data = "ExtraDimensionsLED:ffbar2gammagamma = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "off")
{
$data = "ExtraDimensionsLED:gg2gammagamma = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "off")
{
$data = "ExtraDimensionsLED:Trunc = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "2")
{
$data = "ExtraDimensionsLED:n = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "2000.")
{
$data = "ExtraDimensionsLED:MD = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "2000.")
{
$data = "ExtraDimensionsLED:LambdaT = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "off")
{
$data = "ExtraDimensionsUnpart:ffbar2UZ = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
if($_POST["21"] != "off")
{
$data = "ExtraDimensionsUnpart:ffbar2Ugamma = ".$_POST["21"]."\n";
fwrite($handle,$data);
}
if($_POST["22"] != "off")
{
$data = "ExtraDimensionsUnpart:ffbar2gammagamma = ".$_POST["22"]."\n";
fwrite($handle,$data);
}
if($_POST["23"] != "off")
{
$data = "ExtraDimensionsUnpart:gg2gammagamma = ".$_POST["23"]."\n";
fwrite($handle,$data);
}
if($_POST["24"] != "off")
{
$data = "ExtraDimensionsUnpart:Trunc = ".$_POST["24"]."\n";
fwrite($handle,$data);
}
if($_POST["25"] != "2")
{
$data = "ExtraDimensionsUnpart:spinU = ".$_POST["25"]."\n";
fwrite($handle,$data);
}
if($_POST["26"] != "1.4")
{
$data = "ExtraDimensionsUnpart:dU = ".$_POST["26"]."\n";
fwrite($handle,$data);
}
if($_POST["27"] != "2000.")
{
$data = "ExtraDimensionsUnpart:LambdaU = ".$_POST["27"]."\n";
fwrite($handle,$data);
}
if($_POST["28"] != "1.0")
{
$data = "ExtraDimensionsUnpart:lambda = ".$_POST["28"]."\n";
fwrite($handle,$data);
}
if($_POST["29"] != "1.0")
{
$data = "ExtraDimensionsUnpart:ratio = ".$_POST["29"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2009 Torbjorn Sjostrand -->

