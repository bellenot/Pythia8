<html>
<head>
<title>Electroweak Processes</title>
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

<form method='post' action='ElectroweakProcesses.php'>

<h2>Electroweak Processes</h2>

This page contains processes involving Prompt-photon, <i>gamma^*/Z^0</i> 
and <i>W^+-</i> production, plus a few with <i>t</i>-channel boson 
exchange. Many of the processes are not yet completed, e.g. there are not
yet the correct angular decay distributions in many cases.

<h3>Prompt photon processes</h3>

This group collects the processes where one or two photons are
produced by the hard process. Additional sources of photons 
include parton showers and hadron decays. A <i>pT</i> cut
is required to stay away from the unphysical low-<i>pT</i> region.

<br/><br/><strong>PromptPhoton:all</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of all prompt photon processes, 
as listed separately in the following.
  

<br/><br/><strong>PromptPhoton:qg2qgamma</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q g -> q gamma</i>.
Code 201.
  

<br/><br/><strong>PromptPhoton:qqbar2ggamma</strong>  <input type="radio" name="3" value="on"><strong>On</strong>
<input type="radio" name="3" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar -> g gamma</i>.
Code 202.
  

<br/><br/><strong>PromptPhoton:gg2ggamma</strong>  <input type="radio" name="4" value="on"><strong>On</strong>
<input type="radio" name="4" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> g gamma</i>.
<br/><b>Note:</b> This is a box graph. The full quark-mass 
dependence in the loop leads to very complicated expressions. 
The current implementation is based on assuming five massless 
quarks, and thus is questionable at small (<i>pT &lt; m_b</i>) 
or large (<i>pT > m_t</i>) transverse momenta.
Code 203.
  

<br/><br/><strong>PromptPhoton:qqbar2gammagamma</strong>  <input type="radio" name="5" value="on"><strong>On</strong>
<input type="radio" name="5" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>q qbar -> gamma gamma</i>.
Code 204.
  

<br/><br/><strong>PromptPhoton:gg2gammagamma</strong>  <input type="radio" name="6" value="on"><strong>On</strong>
<input type="radio" name="6" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>g g -> gamma gamma</i>.
<br/><b>Note:</b> This is a box graph. The full quark-mass 
dependence in the loop leads to very complicated expressions. 
The current implementation is based on assuming five massless 
quarks, and thus is questionable at small (<i>pT &lt; m_b</i>) 
or large (<i>pT > m_t</i>) transverse momenta.
Code 205.
  

<h3>Weak boson processes</h3>

Under this heading we group processes involving the production
of a single electroweak gauge boson, i.e. a <i>gamma^*/Z^0</i>
or a <i>W^+-</i>, or a pair of them, or one of them in 
combination with a parton. Since the three sets are partly 
conflicting, each is associated with its own group flag.
In addition, <i>t</i>-channel exchange of such a boson 
between two fermions form a separate group.

<h4>Boson exchange</h4>

<br/><br/><strong>WeakBosonExchange:all</strong>  <input type="radio" name="7" value="on"><strong>On</strong>
<input type="radio" name="7" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of <i>gamma^*/Z^0</i>
or <i>W^+-</i> exchange between two fermions.
  

<br/><br/><strong>WeakBosonExchange:ff2ff(t:gmZ)</strong>  <input type="radio" name="8" value="on"><strong>On</strong>
<input type="radio" name="8" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f f' -> f f'</i> via </i>gamma^*/Z^0</i>
<i>t</i>-channel exchange, with full interference
between the <i>gamma^*</i> and <i>Z^0</i>.
Code 211.
  

<br/><br/><strong>WeakBosonExchange:ff2ff(t:W)</strong>  <input type="radio" name="9" value="on"><strong>On</strong>
<input type="radio" name="9" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f_1 f_2 -> f_3 f_4</i> via </i>W^+-</i>
<i>t</i>-channel exchange.
Code 212.
  

<h4>Single boson</h4>

<br/><br/><strong>WeakSingleBoson:all</strong>  <input type="radio" name="10" value="on"><strong>On</strong>
<input type="radio" name="10" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of a single <i>gamma^*/Z^0</i>
or <i>W^+-</i> production.
  

<br/><br/><strong>WeakSingleBoson:ffbar2gmZ</strong>  <input type="radio" name="11" value="on"><strong>On</strong>
<input type="radio" name="11" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> gamma^*/Z^0</i>, with full interference
between the <i>gamma^*</i> and <i>Z^0</i>.
Code 221.
  

<br/><br/><strong>WeakSingleBoson:ffbar2W</strong>  <input type="radio" name="12" value="on"><strong>On</strong>
<input type="radio" name="12" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar' -> W^+-</i>.
Code 222.
  

<br/><br/><strong>WeakSingleBoson:ffbar2ffbar(s:gm)</strong>  <input type="radio" name="13" value="on"><strong>On</strong>
<input type="radio" name="13" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> gamma^* -> f' fbar'</i>. Subset of 
process 221, but written as a <i>2 -> 2</i> process, so that 
<i>pT</i> can be used as ordering variable, e.g. in multiple 
interactions. As a consequence the scale for pdf's is <i>pT</i>
rather than <i>m</i>, so results will not be identical even 
for the same phase space cuts. Hardcoded for the final state 
being either of the five quark flavours or three lepton ones. 
Not included in the <code>WeakSingleBoson:all</code> set. 
Code 223.
  

<h4>Boson pair</h4>

The processes in this subset are not yet fully handled.

<br/><br/><strong>WeakDoubleBoson:all</strong>  <input type="radio" name="14" value="on"><strong>On</strong>
<input type="radio" name="14" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of pair production of <i>gamma^*/Z^0</i>
and <i>W^+-</i>.
  
 
<br/><br/><strong>WeakDoubleBoson:ffbar2ZW</strong>  <input type="radio" name="15" value="on"><strong>On</strong>
<input type="radio" name="15" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar' -> Z^0 W^+-</i>. Note that here the 
<i>gamma^*</i> is not included.
Code 232.
  
 
<br/><br/><strong>WeakDoubleBoson:ffbar2WW</strong>  <input type="radio" name="16" value="on"><strong>On</strong>
<input type="radio" name="16" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> W^+ W^-</i>.
Code 233.
  

<h4>Boson and parton</h4>

The processes in this subset are not yet fully handled.

<br/><br/><strong>WeakBosonAndParton:all</strong>  <input type="radio" name="17" value="on"><strong>On</strong>
<input type="radio" name="17" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for the group of production of a single electroweak 
gauge boson, i.e. a <i>gamma^*/Z^0</i> or a <i>W^+-</i>, in 
association with a parton, i.e. a quark, gluon, photon or lepton.
These processes give first-order corrections to the ones in the
<code>WeakSingleBoson</code> class, and boths sets cannot be used
simultaneously without unphysical doublecounting. The current class
should only be used to study the high-<i>pT</i> tail of the 
gauge-boson production processes, while the ones in 
<code>WeakSingleBoson</code> should be used for inclusive production.
  
 
<br/><br/><strong>WeakBosonAndParton:qqbar2Wg</strong>  <input type="radio" name="18" value="on"><strong>On</strong>
<input type="radio" name="18" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> W^+- g</i>.
Code 246.
  
 
<br/><br/><strong>WeakBosonAndParton:qg2Wq</strong>  <input type="radio" name="19" value="on"><strong>On</strong>
<input type="radio" name="19" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f g -> f W^+-</i>.
Code 247.
  
 
<br/><br/><strong>WeakBosonAndParton:ffbar2Wgm</strong>  <input type="radio" name="20" value="on"><strong>On</strong>
<input type="radio" name="20" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Scattering <i>f fbar -> W^+- gamma</i>.
Code 248.
  

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
$data = "PromptPhoton:all = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "PromptPhoton:qg2qgamma = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "off")
{
$data = "PromptPhoton:qqbar2ggamma = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "off")
{
$data = "PromptPhoton:gg2ggamma = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "off")
{
$data = "PromptPhoton:qqbar2gammagamma = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "off")
{
$data = "PromptPhoton:gg2gammagamma = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "off")
{
$data = "WeakBosonExchange:all = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "off")
{
$data = "WeakBosonExchange:ff2ff(t:gmZ) = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "off")
{
$data = "WeakBosonExchange:ff2ff(t:W) = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "off")
{
$data = "WeakSingleBoson:all = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "off")
{
$data = "WeakSingleBoson:ffbar2gmZ = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "off")
{
$data = "WeakSingleBoson:ffbar2W = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "off")
{
$data = "WeakSingleBoson:ffbar2ffbar(s:gm) = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "off")
{
$data = "WeakDoubleBoson:all = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "off")
{
$data = "WeakDoubleBoson:ffbar2ZW = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "off")
{
$data = "WeakDoubleBoson:ffbar2WW = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "off")
{
$data = "WeakBosonAndParton:all = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "off")
{
$data = "WeakBosonAndParton:qqbar2Wg = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "off")
{
$data = "WeakBosonAndParton:qg2Wq = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "off")
{
$data = "WeakBosonAndParton:ffbar2Wgm = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright C 2007 Torbjorn Sjostrand -->
