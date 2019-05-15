<html>
<head>
<title>PDF Selection</title>
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

<form method='post' action='PDFSelection.php'>

<h2>PDF Selection</h2>

This page contains three subsections. The first deals with how to 
pick  the parton distribution set for protons, including from LHAPDF, 
to be used for all proton and antiproton beams. The second is a special
option that allows a separate PDF set to be used for the hard process
only, while the first choice would still apply to everything else.
The third gives the possibility to switch off the lepton 
"parton density". 

<h3>Parton densities for protons</h3>

There is one main physics choice to be made with the <code>Pythia</code> 
class, namely which parton densities to use, a choice that then is 
propagated through the program. 
<br/><b>Warning 1:</b> the choice of PDF set affects a number of
properties of events. A change of PDF therefore requires a complete 
retuning e.g.  of the multiple-interactions model for minimum-bias and 
underlying events.
<br/><b>Warning 2:</b> People often underestimate the differences 
between different sets on the market. The sets are constructed to behave 
more or less similarly at large <i>x</i> and <i>Q2</i>, while the 
multiple interactions are dominated by the behaviour in the region of 
small <i>x</i> and <i>Q2</i>. A good PDF parametrization ought to be
sensible down to <i>x = 10^{-6}</i> (<i>x = 10^{-7}</i>) and
<i>Q2 = 1</i> GeV2 for Tevatron (LHC) applications. Unfortunately there
are distributions on the market that completely derail in that region.
The <code>main41.cc</code> and <code>main42.cc</code> programs in the 
<code>examples</code> subdirectory provide some examples of absolutely
minimal sanity checks before a new PDF set is put in production.
<br/><b>Warning 3:</b> Do not blindly assume that an NLO tune has to be 
better than an LO one when combined with the LO matrix elements in PYTHIA.
There are explicit examples where such thinking can lead you down the 
wrong alley.

<p/>
The simplest option is to pick one 
of the few distributions available internally:

<br/><br/><table><tr><td><strong>PDF:pSet  </td><td>  &nbsp;&nbsp;(<code>default = <strong>2</strong></code>; <code>minimum = 1</code>; <code>maximum = 2</code>)</td></tr></table>
<modepick name="PDF:pSet" default="2" min="1" max="2">
Parton densities to be used for proton beams (and, by implication,
antiproton ones):
<br/>
<input type="radio" name="1" value="1"><strong>1 </strong>: GRV 94 L;<br/>
<input type="radio" name="1" value="2" checked="checked"><strong>2 </strong>: CTEQ 5 L.<br/>
</modepick> 

<p/>
Obviously this choice is mainly intended to get going, and if you link to
the <a href="http://projects.hepforge.org/lhapdf/" target="page">LHAPDF 
library</a> [<a href="Bibliography.php" target="page">Wha05</a>] you get access to a much wider selection.

<br/><br/><strong>PDF:useLHAPDF</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
If off then the choice of proton PDF is based on <code>pPDFset</code>
above. If on then it is instead based on the choice of 
<code>LHAPDFset</code> and <code>LHAPDFmember</code> below.
<br/><b>Note:</b> in order for this option to work you must have 
compiled PYTHIA appropriately and have set the <code>LHAPATH</code> 
environment variable to provide the data-files directory of your local 
LHAPDF installation. See the README file in the <code>examples</code> 
directory for further instructions. 
  

<br/><br/><table><tr><td><strong>PDF:LHAPDFset  </td><td></td><td> <input type="text" name="3" value="MRST2004FF4lo.LHgrid" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>MRST2004FF4lo.LHgrid</strong></code>)</td></tr></table>
Name of proton PDF set from LHAPDF to be used. You have to choose 
from the 
<a href="http://projects.hepforge.org/lhapdf/pdfsets" target="page">
list of available sets</a>. Examples of some recent ones would be 
cteq61.LHpdf, cteq61.LHgrid, cteq6l.LHpdf, cteq6ll.LHpdf, 
MRST2004nlo.LHpdf, MRST2004nlo.LHgrid, MRST2004nnlo.LHgrid and 
MRST2004FF3lo.LHgrid. If you pick a LHpdf set it will require some 
calculation the first time it is called. 
<br/><b>Technical note:</b> if you provide a name beginning with a 
slash (/) it is assumed you want to provide the full file path and then
<code>initPDFsetM(name)</code> is called, else the correct path is assumed 
already set and <code>initPDFsetByNameM(name)</code> is called.
   

<br/><br/><table><tr><td><strong>PDF:LHAPDFmember  </td><td></td><td> <input type="text" name="4" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>)</td></tr></table>
Further choice of a specific member from the set picked above. Member 0
should normally correspond to the central value, with higher values
corresponding to different error PDF's somewhat off in different 
directions. You have to check from set to set which options are open.
<br/><b>Note:</b> you can only use one member in a run, so if you
want to sweep over many members you either have to do many separate
runs or, as a simplification, save the 
<?php $filepath = $_GET["filepath"];
echo "<a href='EventInformation.php?filepath=".$filepath."' target='page'>";?>pdf weights</a> at the hard scattering
and do an offline reweighting of events.
</modeopen>   

<p/> 
The current LHAPDF code does not provide a way to find the <i>x</i>
and <i>Q2</i> limits of a set, nor does it guarantee a sensible 
behaviour outside of those limits. (The LHAGLUE interface fixes this
problem, but is not used here. The limits are reported in the
<a href="http://projects.hepforge.org/lhapdf/manual#tth_sEcA" target="page">
PDF set list</a>, however.) Error messages may abound, or execution 
may stop unexpectedly. In a near future it will become possible to 
interrogate the limits. Meanwhile you can yourself require that 
PDF's should only be invoked in a specified range, to avoid problems.

<br/><br/><strong>PDF:limitLHAPDF</strong>  <input type="radio" name="5" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="5" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
If on then you can set <i>x</i> and <i>Q2</i> limits below.
  

<br/><br/><table><tr><td><strong>PDF:xMinLHAPDF </td><td></td><td> <input type="text" name="6" value="1e-6" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1e-6</strong></code>)</td></tr></table>
Freeze <i>x f_i(x, Q2)</i> below this <i>x</i> value.
  

<br/><br/><table><tr><td><strong>PDF:xMaxLHAPDF </td><td></td><td> <input type="text" name="7" value="0.9999" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.9999</strong></code>)</td></tr></table>
Set <i>x f_i(x, Q2) = 0</i> above this <i>x</i> value.
  

<br/><br/><table><tr><td><strong>PDF:Q2MinLHAPDF </td><td></td><td> <input type="text" name="8" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
Freeze <i>x f_i(x, Q2)</i> below this <i>Q2</i> value.
  

<br/><br/><table><tr><td><strong>PDF:Q2MaxLHAPDF </td><td></td><td> <input type="text" name="9" value="1e8" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1e8</strong></code>)</td></tr></table>
Freeze <i>x f_i(x, Q2)</i> above this <i>Q2</i> value.
  

<p/> 
If you want to use PDF's not found in LHAPDF, or you want to interface
LHAPDF another way, you have full freedom to use the more generic 
<?php $filepath = $_GET["filepath"];
echo "<a href='PartonDistributions.php?filepath=".$filepath."' target='page'>";?>interface options</a>.

<h3>Parton densities for protons in the hard process</h3>

The above options provides a PDF set that will be used everywhere:
for the hard process, the parton showers and the multiple interactions
alike. As already mentioned, therefore a change of PDF should be
accompanied by a <b>complete</b> retuning of the whole MI framework,
and maybe more. There are cases where one may want to explore 
different PDF options for the hard process, but would not want to touch 
the rest. If several different sets are to be compared, a simple
reweighting based on the <?php $filepath = $_GET["filepath"];
echo "<a href='EventInformation.php?filepath=".$filepath."' target='page'>";?>originally 
used</a> flavour, <i>x</i>, <i>Q2</i> and PDF values may offer the 
best route. The options in this section allow a choice of the PDF set
for the hard process alone, while the choice made in the previous section
would still be used for everything else. The hardest interaction
of the minimum-bias process is part of the multiple-interactions
framework and so does not count as a hard process here. 

<p/>
Of course it is inconsistent to use different PDF's in different parts 
of an event, but if the <i>x</i> and <i>Q2</i> ranges mainly accessed 
by the components are rather different then the contradiction would not be
too glaring. Furthermore, since standard PDF's are one-particle-inclusive
we anyway have to 'invent' our own PDF modifications to handle configurations
where more than one parton is kicked out of the proton [<a href="Bibliography.php" target="page">Sjo04</a>]. 

<p/>
The PDF choices that can be made are the same as above, so we do not 
repeat the detailed discussion.

<br/><br/><strong>PDF:useHard</strong>  <input type="radio" name="10" value="on"><strong>On</strong>
<input type="radio" name="10" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
If on then select a separate PDF set for the hard process, using the 
variables below. If off then use the same PDF set for everything,
as already chosen above.   
  

<br/><br/><table><tr><td><strong>PDF:pHardSet  </td><td>  &nbsp;&nbsp;(<code>default = <strong>2</strong></code>; <code>minimum = 1</code>; <code>maximum = 2</code>)</td></tr></table>
<modepick name="PDF:pHardSet" default="2" min="1" max="2">
Parton densities to be used for proton beams (and, by implication,
antiproton ones):
<br/>
<input type="radio" name="11" value="1"><strong>1 </strong>: GRV 94 L;<br/>
<input type="radio" name="11" value="2" checked="checked"><strong>2 </strong>: CTEQ 5 L.<br/>
</modepick> 

<br/><br/><strong>PDF:useHardLHAPDF</strong>  <input type="radio" name="12" value="on"><strong>On</strong>
<input type="radio" name="12" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
If off then the choice of proton PDF is based on <code>hardpPDFset</code>
above. If on then it is instead based on the choice of 
<code>hardLHAPDFset</code> and <code>hardLHAPDFmember</code> below.
  

<br/><br/><table><tr><td><strong>PDF:hardLHAPDFset  </td><td></td><td> <input type="text" name="13" value="MRST2004FF4lo.LHgrid" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>MRST2004FF4lo.LHgrid</strong></code>)</td></tr></table>
Name of proton PDF set from LHAPDF to be used. 
   

<br/><br/><table><tr><td><strong>PDF:hardLHAPDFmember  </td><td></td><td> <input type="text" name="14" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>)</td></tr></table>
Further choice of a specific member from the set picked above. 
</modeopen>   

<br/><br/><strong>PDF:limitHardLHAPDF</strong>  <input type="radio" name="15" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="15" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
If on then you can set <i>x</i> and <i>Q2</i> limits below.
  

<br/><br/><table><tr><td><strong>PDF:xMinHardLHAPDF </td><td></td><td> <input type="text" name="16" value="1e-6" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1e-6</strong></code>)</td></tr></table>
Freeze <i>x f_i(x, Q2)</i> below this <i>x</i> value.
  

<br/><br/><table><tr><td><strong>PDF:xMaxHardLHAPDF </td><td></td><td> <input type="text" name="17" value="0.9999" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.9999</strong></code>)</td></tr></table>
Set <i>x f_i(x, Q2) = 0</i> above this <i>x</i> value.
  

<br/><br/><table><tr><td><strong>PDF:Q2MinHardLHAPDF </td><td></td><td> <input type="text" name="18" value="1." size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.</strong></code>)</td></tr></table>
Freeze <i>x f_i(x, Q2)</i> below this <i>Q2</i> value.
  

<br/><br/><table><tr><td><strong>PDF:Q2MaxHardLHAPDF </td><td></td><td> <input type="text" name="19" value="1e8" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1e8</strong></code>)</td></tr></table>
Freeze <i>x f_i(x, Q2)</i> above this <i>Q2</i> value.
  

<h3>Parton densities for leptons</h3>

For electrons/leptons there is no need to choose between different 
parametrizations, since only one implementation is available, and 
should be rather uncontroversial (apart from some technical details).
However, insofar as e.g. <i>e^+ e^-</i> data often are corrected 
back to a world without any initial-state photon radiation, it is 
useful to have a corresponding option available here.

<br/><br/><strong>PDF:lepton</strong>  <input type="radio" name="20" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="20" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
Use parton densities for lepton beams or not. If off the colliding
leptons carry the full beam energy, if on part of the energy is 
radiated away by initial-state photons. In the latter case the
initial-state showers will generate the angles and energies of the
set of photons that go with the collision. In addition one collinear
photon per beam carries any leftover amount of energy not described
by shower emissions. If the initial-state showers are switched off 
these collinear photons will carry the full radiated energy.  
   

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

if($_POST["1"] != "2")
{
$data = "PDF:pSet = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "PDF:useLHAPDF = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "MRST2004FF4lo.LHgrid")
{
$data = "PDF:LHAPDFset = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "0")
{
$data = "PDF:LHAPDFmember = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "on")
{
$data = "PDF:limitLHAPDF = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "1e-6")
{
$data = "PDF:xMinLHAPDF = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "0.9999")
{
$data = "PDF:xMaxLHAPDF = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "1.")
{
$data = "PDF:Q2MinLHAPDF = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "1e8")
{
$data = "PDF:Q2MaxLHAPDF = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "off")
{
$data = "PDF:useHard = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "2")
{
$data = "PDF:pHardSet = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "off")
{
$data = "PDF:useHardLHAPDF = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "MRST2004FF4lo.LHgrid")
{
$data = "PDF:hardLHAPDFset = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "0")
{
$data = "PDF:hardLHAPDFmember = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "on")
{
$data = "PDF:limitHardLHAPDF = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "1e-6")
{
$data = "PDF:xMinHardLHAPDF = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "0.9999")
{
$data = "PDF:xMaxHardLHAPDF = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "1.")
{
$data = "PDF:Q2MinHardLHAPDF = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
if($_POST["19"] != "1e8")
{
$data = "PDF:Q2MaxHardLHAPDF = ".$_POST["19"]."\n";
fwrite($handle,$data);
}
if($_POST["20"] != "on")
{
$data = "PDF:lepton = ".$_POST["20"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2007 Torbjorn Sjostrand -->
