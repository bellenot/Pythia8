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

This page contains Higgs production in the Standard Model.
It should eventually be expanded at least also to cover the MSSM. 

<h3>Standard-Model Higgs, basic processes</h3>

This section provides the standard set of processes that can be
run together to provide a reasonably complete overview of possible
production channels for a single Standard-Model Higgs.

<br/><br/><strong>SMHiggs:all</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of Higgs production within the Standard Model.
  

<br/><br/><strong>SMHiggs:ffbar2H</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> H^0</i>, where <i>f</i> sums over available
flavours except top. Related to the mass-dependent Higgs point coupling 
to fermions, so at hadron colliders the bottom contribution will
dominate.
Code 901.
  

<br/><br/><strong>SMHiggs:gg2H</strong>  <input type="radio" name="3" value="on"><strong>On</strong>
<input type="radio" name="3" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> H^0</i> via loop contributions primarily from
top.
Code 902.
  

<br/><br/><strong>SMHiggs:gmgm2H</strong>  <input type="radio" name="4" value="on"><strong>On</strong>
<input type="radio" name="4" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>gamma gamma -> H^0</i> via loop contributions primarily 
from top and <i>W</i>.
Code 903.
  

<br/><br/><strong>SMHiggs:ffbar2HZ</strong>  <input type="radio" name="5" value="on"><strong>On</strong>
<input type="radio" name="5" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> H^0 Z^0</i> via <i>s</i>-channel <i>Z^0</i>
exchange.
Code 904.
  

<br/><br/><strong>SMHiggs:ffbar2HW</strong>  <input type="radio" name="6" value="on"><strong>On</strong>
<input type="radio" name="6" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> H^0 W^+-</i> via <i>s</i>-channel <i>W^+-</i>
exchange.
Code 905.
  

<br/><br/><strong>SMHiggs:ff2Hff(t:ZZ)</strong>  <input type="radio" name="7" value="on"><strong>On</strong>
<input type="radio" name="7" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f f' -> H^0 f f'</i> via <i>Z^0 Z^0</i> fusion.
Code 906.
  

<br/><br/><strong>SMHiggs:ff2Hff(t:WW)</strong>  <input type="radio" name="8" value="on"><strong>On</strong>
<input type="radio" name="8" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f_1 f_2 -> H^0 f_3 f_4</i> via <i>W^+ W^-</i> fusion.
Code 907.
  

<br/><br/><strong>SMHiggs:gg2Httbar</strong>  <input type="radio" name="9" value="on"><strong>On</strong>
<input type="radio" name="9" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> H^0 t tbar</i> via <i>t tbar</i> fusion
(or, alternatively put, Higgs radiation off a top line).
Warning: unfortunately this process is rather slow, owing to a
lengthy cross-section expression and inefficient phase-space selection.
Code 908.
  

<br/><br/><strong>SMHiggs:qqbar2Httbar</strong>  <input type="radio" name="10" value="on"><strong>On</strong>
<input type="radio" name="10" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar -> H^0 t tbar</i> via <i>t tbar</i> fusion
(or, alternatively put, Higgs radiation off a top line).
Warning: unfortunately this process is rather slow, owing to a
lengthy cross-section expression and inefficient phase-space selection.
Code 909.
  

<h3>Standard-Model Higgs, further processes</h3>

A number of further production processes has been implemented, that 
are specializations of some of the above ones to the high-<i>pT</i> 
region. The sets therefore could not be used simultaneously
without unphysical doublecounting, as further explained below. 
They are not switched on by the <code>SMHiggs:all</code> flag, but 
have to be switched on for each separate process after due consideration.

<p/>
The first three processes in this section are related to the Higgs
point coupling to fermions, and so primarily are of interest for 
<i>b</i> quarks. It is here useful to begin by reminding that 
a process like <i>b bbar -> H^0</i> implies that a <i>b/bbar</i> 
is taken from each incoming hadron, leaving behind its respective
antiparticle. The initial-state showers will then add one 
<i>g -> b bbar</i> branching on either side, so that effectively
the process becomes <i>g g -> H0 b bbar</i>. This would be the
same basic process as the <i>g g -> H^0 t tbar</i> one used for top.
The difference is that (a) no PDF's are defined for top and 
(b) the shower approach would not be good enough to provide sensible
kinematics for the <i>H^0 t tbar</i> subsystem. By contrast, owing 
to the <i>b</i> being much lighter than the Higgs, multiple 
gluon emissions must be resummed for <i>b</i>, as is done by PDF's 
and showers, in order to obtain a sensible description of the total 
production rate,  when the <i>b</i> quarks predominantly are produced 
at small <i>pT</i> values.

<br/><br/><strong>SMHiggs:qg2Hq</strong>  <input type="radio" name="11" value="on"><strong>On</strong>
<input type="radio" name="11" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>b g -> H^0 b</i>. This process gives first-order 
corrections to the <i>f fbar -> H^0</i> one above, and should only be 
used to study  the high-<i>pT</i> tail, while <i>f fbar -> H^0</i> 
should be used for inclusive production. Only the dominant <i>c</i>c 
and <i>b</i> contributions are included, and generated separately 
for technical reasons. 
Code 911.

  
<br/><br/><strong>SMHiggs:gg2Hbbbar</strong>  <input type="radio" name="12" value="on"><strong>On</strong>
<input type="radio" name="12" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> H^0 b bbar</i>. This process is yet one order 
higher of the <i>b bbar -> H^0</i> and <i>b g -> H^0 b</i> chain,
where now two quarks should be required above some large <i>pT</i>
threshold.
Warning: unfortunately this process is rather slow, owing to a
lengthy cross-section expression and inefficient phase-space selection.
Code 912.
  

  
<br/><br/><strong>SMHiggs:qqbar2Hbbbar</strong>  <input type="radio" name="13" value="on"><strong>On</strong>
<input type="radio" name="13" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar -> H^0 b bbar</i> via an <i>s</i>-channel
gluon, so closely related to the previous one, but typically less 
important owing to the smaller rate of (anti)quarks relative to 
gluons.
Warning: unfortunately this process is rather slow, owing to a
lengthy cross-section expression and inefficient phase-space selection.
Code 913.
  

<p/>
The second set of processes are predominantly first-order corrections 
to the <i>g g -> H^0</i> process, again dominated by the top loop.
We here only provide the kinematical expressions obtained in the 
limit that the top quark goes to infinity, but scaled to the 
finite-top-mass coupling in <i>g g -> H^0</i>. (Complete loop
expressions are available e.g. in PYTHIA 6.4 but are very lengthy.) 
This provides a reasonably accurate description for "intermediate" 
<i>pT</i> values, but fails when the <i>pT</i> scale approaches
the top mass. 
 
<br/><br/><strong>SMHiggs:gg2Hg(l:t)</strong>  <input type="radio" name="14" value="on"><strong>On</strong>
<input type="radio" name="14" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> H^0 g</i> via loop contributions primarily 
from top.
Code 914.
  
 
<br/><br/><strong>SMHiggs:qg2Hq(l:t)</strong>  <input type="radio" name="15" value="on"><strong>On</strong>
<input type="radio" name="15" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q g -> H^0 q</i> via loop contributions primarily 
from top. Not to be confused with the <code>SMHiggs:bg2Hb</code>
process above, with its direct fermion-to-Higgs coupling.
Code 915.
  
 
<br/><br/><strong>SMHiggs:qqbar2Hg(l:t)</strong>  <input type="radio" name="16" value="on"><strong>On</strong>
<input type="radio" name="16" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar -> H^0 g</i> via an <i>s</i>-channel gluon
and loop contributions primarily from top. Is strictly speaking a 
"new" process, not directly derived from <i>g g -> H^0</i>, and
could therefore be included in the standard mix without doublecounting, 
but is numerically negligible.
Code 916.
  


<h3>Parameters for Higgs production and decay</h3>

<br/><br/><strong>ResonanceSMH:linearWidthWWZZ</strong>  <input type="radio" name="17" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="17" value="off"><strong>Off</strong>
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
  

<br/><br/><table><tr><td><strong>ResonanceSMH:parity  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>; <code>maximum = 3</code>)</td></tr></table>
<modepick name="ResonanceSMH:parity" default="1" min="0" max="3">
possibility to modify angular decay correlations in Higgs decay to
<ei>Z^0 Z^0</ei> or <ei>W^+ W^-</ei> to four fermions. Currently it 
does not affect the partial width of the channels.
<br/>
<input type="radio" name="18" value="0"><strong>0 </strong>: isotropic decays.<br/>
<input type="radio" name="18" value="1" checked="checked"><strong>1 </strong>: assuming the Higgs is a pure scalar (CP-even), as in the Standard Model.<br/>
<input type="radio" name="18" value="2"><strong>2 </strong>: assuming the Higgs is a pure pseudoscalar(CP-odd).<br/>
<input type="radio" name="18" value="3"><strong>3 </strong>: assuming the Higgs is a mixture of the two, including the CP-violating interference term. The parameter<ei>eta</ei>, see below, sets the strength of the CP-odd admixture,with the interference term being proportional to <ei>eta</ei>and the CP-odd one to <ei>eta^2</ei>.<br/>
</modepick>

<br/><br/><table><tr><td><strong>ResonanceSMH:etaParity </td><td></td><td> <input type="text" name="19" value="0." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.</strong></code>)</td></tr></table>
The <i>eta</i> value of CP-violation in the 
<code>ResonanceSMH:parity = 2</code> option. 
  

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
$data = "SMHiggs:ff2Hff(t:ZZ) = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "off")
{
$data = "SMHiggs:ff2Hff(t:WW) = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "off")
{
$data = "SMHiggs:gg2Httbar = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "off")
{
$data = "SMHiggs:qqbar2Httbar = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "off")
{
$data = "SMHiggs:qg2Hq = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "off")
{
$data = "SMHiggs:gg2Hbbbar = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "off")
{
$data = "SMHiggs:qqbar2Hbbbar = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "off")
{
$data = "SMHiggs:gg2Hg(l:t) = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "off")
{
$data = "SMHiggs:qg2Hq(l:t) = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "off")
{
$data = "SMHiggs:qqbar2Hg(l:t) = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "on")
{
$data = "ResonanceSMH:linearWidthWWZZ = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "1")
{
$data = "ResonanceSMH:parity = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "0.")
{
$data = "ResonanceSMH:etaParity = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2007 Torbjorn Sjostrand -->

