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

The implementation of SUSY processes is currently still in a <b>development</b>
stage, so careful case-by-case validations against other codes are
strongly recommended and should be considered mandatory. 
Please notify the authors (skands@fnal.gov) 
of any significant deviations (i.e., larger than 10%). In some cases,
limited preliminary validations have 
already been carried out by the authors. This is remarked on process by
process below. 

<p/>
Here is collected processes involving supersymmetric particle 
production, with the exception of the (extended) Higgs sector.
Since the number of separate but closely related processes is so big,
there will not be switches for each separate process but only for a
reasonable set of subgroups. However, the general 
switches <code>SUSY:idA</code> and <code>SUSY:idB</code> may be used in
conjunction with any of these groups to provide some additional
flexibility to concentrate on processes involving only specific (s)particle
final states, see below. 

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

<p/>
Note also that lepton and photon initial states are not yet available. Only
quark/gluon-initiated 2->2 processes have been implemented. Likewise, direct
slepton production has not yet been implemented (i.e., 2->2 processes
involving sleptons in the final state). Sleptons will of course still be
produced through cascade decays of heavier (s)particles. 

<p/>
Finally, note that these cross sections will be correctly folded with open
branching fractions of cascade decays, but at present any difference between
particle and antiparticle decay tables is not taken into account. This
possibility will be included in a future update. 

<br/><br/><strong>SUSY:all</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Common switch for production of supersymmetric particles, i.e. 
particles with R-parity -1. 
  

<br/><br/><table><tr><td><strong>SUSY:idA  </td><td></td><td> <input type="text" name="2" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>)</td></tr></table>
Option to limit the sum over possible outgoing states in SUSY
<i>2 -> 2</i> processes to ones including a specific particle 
identity code. The default corresponds to summing over all possible 
indices. A non-zero value of <code>SUSY:idA</code> selects only processes 
that contain the state corresponding to that particular particle identity 
code in the fundamental <i>2 -> 2</i> scattering process (symmetrized 
over particle/antiparticle). It is the user's responsibility to ensure 
that (a subset of) the processes be to simulated actually include this 
particle at the <i>2 -> 2</i> level; thus, asking for the lightest 
neutralino (code 1000021) to be present in a squark-squark production 
process will give no match. 
  

<br/><br/><table><tr><td><strong>SUSY:idB  </td><td></td><td> <input type="text" name="3" value="0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>)</td></tr></table>
As for <code>SUSY:idA</code>, but requires an additional particle
with PDG code <code>SUSY:idB</code> to be present in the 2->2
process. Thus, using  <code>SUSY:idA</code> and  <code>SUSY:idB</code> 
a specific subprocess can be selected. Again only the absolute sign is 
used, i.e. the summation over particle and antiparticle is retained.
  

<h3>Gluino Pair Production</h3>

<br/><br/><strong>SUSY:gg2gluinogluino</strong>  <input type="radio" name="4" value="on"><strong>On</strong>
<input type="radio" name="4" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production of gluinos by gluon-gluon initial states. Does not
depend on flavour violation. This cross section has been validated against
Pythia 6 for SPS1a. 
  

<br/><br/><strong>SUSY:qqbar2gluinogluino</strong>  <input type="radio" name="5" value="on"><strong>On</strong>
<input type="radio" name="5" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production of gluinos by quark-antiquark annihilation and
t-channel squark exchange. So far, these cross sections assume the
squark CKM is aligned with the quark CKM and that all quantities are
real, so effects of non-minimal flavour violation and/or CP violation
are not yet included. This cross section has been validated against
Pythia 6 for SPS1a. 
  

<h3>Associated Squark-Gluino Production</h3>

<br/><br/><strong>SUSY:qg2squarkgluino</strong>  <input type="radio" name="6" value="on"><strong>On</strong>
<input type="radio" name="6" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Associated production of a squark with a gluino. Only implemented for the
flavor-diagonal case. These cross sections have been validated against
Pythia 6 for SPS1a. Note: these cross sections are so far 
limited to the SLHA1 flavor structure. They will not give correct results if
used with an SLHA2 spectrum. 
  

<h3>Squark Pair Production</h3>

<br/><br/><strong>SUSY:gg2squarkantisquark</strong>  <input type="radio" name="7" value="on"><strong>On</strong>
<input type="radio" name="7" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production of a scalar quark together with a scalar antiquark by
gluon annhilation via s-channel gluon exhange, t- and u-channel squark
exchange, and the direct 4-point coupling. 
The cross section expression follows 
[<a href="Bibliography.php" target="page">Boz07</a>] and has been validated against Pythia 6.
  

<br/><br/><strong>SUSY:qqbar2squarkantisquark</strong>  <input type="radio" name="8" value="on"><strong>On</strong>
<input type="radio" name="8" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production of a scalar quark together with a scalar antiquark 
by quark-antiquark annihilation. 
For same-isospin ~q~q* production (i.e., ~u~u*, ~u~c*, ...), 
the s-channel gluon, photon, and Z and
t-channel gluino contributions have so far been implemented (i.e., the
t-channel neutralino contributions are neglected). 
For opposite-isospin ~q~q* production (~u~d*, ~u~s*, ...), 
the s-channel W and t-channel gluino
contributions have been implemented (i.e., the t-channel neutralino
contributions are neglected). 
The cross section expressions follow [<a href="Bibliography.php" target="page">Boz07</a>] 
and should thus be valid also in the case of non-minimal flavour 
violation and/or CP violation. These cross sections have been validated
against Pythia 6 for SPS1a, 
for ~t1~t1*, ~t2~t2*, ~b1~b1*, ~b2~b2*, and ~b1~b2* (+c.c.) production. 
For ~t1~t2* (+c.c.), 
there is currently a factor 2 between Pythia 6 and this implementation; 
pending resolution. The FLV and CPV cases have not yet been validated. 
  

<br/><br/><strong>SUSY:qqbar2squarkantisquark:gluonOnly</strong>  <input type="radio" name="9" value="on"><strong>On</strong>
<input type="radio" name="9" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
When switched <code>on</code> this flag switches off all but the s-channel
gluon contribution in the
calculation of same-isospin squark-antisquark production cross
sections. Intended for reference only. For the most
accurate physics simulation, leave this flag in the <code>off</code>
position.  
  

<br/><br/><strong>SUSY:qq2squarksquark</strong>  <input type="radio" name="10" value="on"><strong>On</strong>
<input type="radio" name="10" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production of scalar quarks (squark-squark and its charge
conjugate process; for squark-antisquark production see above) 
by t- and u-channel gluino, neutralino, and 
 chargino exchange. The cross section expressions follow [<a href="Bibliography.php" target="page">Boz07</a>] 
and should thus be valid also in the case of non-minimal flavour 
violation and/or CP violation. Note: Pythia 6 only included the gluino exchange
 contribution, which typically dominates due to the size of the strong
coupling; for counterchecks, 
the flag <code>SUSY:qq2squarksquark:gluinoOnly</code>
below can be switched on to eliminate the chargino and neutralino
contributions. These cross sections have been validated against Pythia 6 for
SPS1a, for ~b1~b1, ~b1~b2, and ~b2~b2 production, with 5-10% differences
observed.  
  

<br/><br/><strong>SUSY:qq2squarksquark:gluinoOnly</strong>  <input type="radio" name="11" value="on"><strong>On</strong>
<input type="radio" name="11" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
When switched <code>on</code> this flag causes the t- or u-channel 
neutralino and chargino contributions to be ignored in the
calculation of squark pair production cross sections. Intended for reference
only. For the most
accurate physics simulation, leave this flag in the <code>off</code>
position.  
  

<h3>Neutralino and Chargino Pair Production</h3>

<br/><br/><strong>SUSY:qqbar2chi0chi0</strong>  <input type="radio" name="12" value="on"><strong>On</strong>
<input type="radio" name="12" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production of neutralinos by quark-antiquark annihilation. With
four neutralino species this gives ten separate processes, codes 
1201 - 1210. The cross section expressions follow [<a href="Bibliography.php" target="page">Boz07</a>] and
should thus be  valid also in the case of non-minimal flavour
violation and/or CP violation. These cross sections have been validated
against Pythia 6 for SPS1a and against the code 
XSUSY (based on [<a href="Bibliography.php" target="page">Boz07</a>]) for a non-minimal FLV case. 
Validation of CPV cases has not yet been carried out.
  

<br/><br/><strong>SUSY:qqbar2chi+-chi0</strong>  <input type="radio" name="13" value="on"><strong>On</strong>
<input type="radio" name="13" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Associated chargino-neutralino production by quark-antiquark
annihilation. With four neutralino species, two chargino ones, and
maintaining charge conjugate proceeses separate, this gives 16 
separate processes, codes 1221 - 1236. The cross section expressions 
follow [<a href="Bibliography.php" target="page">Boz07</a>] and should thus be valid also in the case of 
non-minimal flavour violation and/or CP violation. These cross sections have
been validated against Pythia 6 for SPS1a and against the code 
XSUSY (based on [<a href="Bibliography.php" target="page">Boz07</a>]) for a non-minimal FLV case. Validation of
 CPV cases has not yet been carried out. 
  

<br/><br/><strong>SUSY:qqbar2chi+chi-</strong>  <input type="radio" name="14" value="on"><strong>On</strong>
<input type="radio" name="14" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production of charginos by quark-antiquark annihilation. With
two chargino species and maintaining mutually charge conjugate
processes separate, this gives four separate processes, codes 
1241 - 1244. The cross section expressions follow [<a href="Bibliography.php" target="page">Boz07</a>] 
and should thus be valid also in the case of non-minimal flavour 
violation and/or CP violation. These cross sections have
been validated 
against Pythia 6 for SPS1a and against the code 
XSUSY (based on [<a href="Bibliography.php" target="page">Boz07</a>]) for a non-minimal FLV case. Validation of
 CPV cases has not yet been carried out. 
  

<h3>Associated Neutralino/Chargino + Squark/Gluino Production</h3>

<br/><br/><strong>SUSY:qg2chi0squark</strong>  <input type="radio" name="15" value="on"><strong>On</strong>
<input type="radio" name="15" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Pair production of neutralinos from quark-gluon initial states.
The cross section expressions follow [<a href="Bibliography.php" target="page">Boz07</a>] and
should thus be  valid also in the case of non-minimal flavour
violation and/or CP violation. Note: these cross sections are in a
development stage and have not yet been validated.
  

<br/><br/><strong>SUSY:qg2chi+-squark</strong>  <input type="radio" name="16" value="on"><strong>On</strong>
<input type="radio" name="16" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Associated chargino-squark production from quark-gluon initial states.
annihilation. The cross section expressions 
follow [<a href="Bibliography.php" target="page">Boz07</a>] and should thus be valid also in the case of 
non-minimal flavour violation and/or CP violation. Note: these cross sections
are in a development stage and have not yet been validated. 
  

<br/><br/><strong>SUSY:qqbar2chi0gluino</strong>  <input type="radio" name="17" value="on"><strong>On</strong>
<input type="radio" name="17" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Associated neutralino-gluino production by quark-antiquark
annihilation. Status: not implemented yet.
  

<br/><br/><strong>SUSY:qqbar2chi+-gluino</strong>  <input type="radio" name="18" value="on"><strong>On</strong>
<input type="radio" name="18" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Associated chargino-gluino production by quark-antiquark
annihilation. Status: not implemented yet.
  

<h3>Slepton Production</h3>

No 2->2 slepton pair production or associated slepton production 
cross sections have been implemented yet. 

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
if($_POST["2"] != "0")
{
$data = "SUSY:idA = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "0")
{
$data = "SUSY:idB = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "off")
{
$data = "SUSY:gg2gluinogluino = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "off")
{
$data = "SUSY:qqbar2gluinogluino = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "off")
{
$data = "SUSY:qg2squarkgluino = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "off")
{
$data = "SUSY:gg2squarkantisquark = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "off")
{
$data = "SUSY:qqbar2squarkantisquark = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "off")
{
$data = "SUSY:qqbar2squarkantisquark:gluonOnly = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "off")
{
$data = "SUSY:qq2squarksquark = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "off")
{
$data = "SUSY:qq2squarksquark:gluinoOnly = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "off")
{
$data = "SUSY:qqbar2chi0chi0 = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "off")
{
$data = "SUSY:qqbar2chi+-chi0 = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "off")
{
$data = "SUSY:qqbar2chi+chi- = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "off")
{
$data = "SUSY:qg2chi0squark = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "off")
{
$data = "SUSY:qg2chi+-squark = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
if($_POST["17"] != "off")
{
$data = "SUSY:qqbar2chi0gluino = ".$_POST["17"]."\n";
fwrite($handle,$data);
}
if($_POST["18"] != "off")
{
$data = "SUSY:qqbar2chi+-gluino = ".$_POST["18"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2009 Peter Skands, Torbjorn Sjostrand -->

