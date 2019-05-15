<html>
<head>
<title>Resonance Decays</title>
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

<form method='post' action='ResonanceDecays.php'>

<h2>Resonance Decays</h2>

The <code>ResonanceDecays</code> class performs the sequential decays of 
all resonances formed in the hard process. Note the important distinction
between "resonances" and other "particles" made in PYTHIA.
<ul> 
<li>
The list of resonances contains <i>gamma^*/Z^0</i>, <i>W^+-</i>, top, 
the Higgs, and essentially all new particles of Beyond-the-Standard-Model 
physics: further Higgses, sfermions, gauginos, techniparticles, and so on. 
The partial widths to different decay channels are perturbatively
calculable, given the parameters of the respective model, and branching
ratios may be allowed to vary across a (reasonably broad) resonance peak.
Usually resonances are short-lived, and therefore it makes sense to consider 
their decays immediately after the primary hard process has been set up. 
Furthermore, in several cases the decay angular distributions are encoded 
as part of the specific process, e.g. the <i>W</i> decays differently in 
<i>f fbar -> W^+-</i>, <i>f fbar -> W^+ W^-</i> and 
<i>h^0 -> W^+ W^- </i>. All of these particles are (in PYTHIA) only 
produced as part of the hard process itself, i.e. they are not produced 
in showers or hadronization processes. Therefore the restriction to 
specific decay channels can be consistently taken into account as a 
corresponding reduction in the cross section of a process. Finally, note 
that all of these resonances have an on-shell mass above 20 GeV, with the 
exception of some hypothetical weakly interacting and stable particles 
such as the gravitino.
</li>
<li>
The other particles include normal hadrons and the Standard-Model leptons,
including the <i>tau^+-</i>. These can be produced in the normal
hadronization and decay description, which involve unknown nonperturbative
parameters and multistep chains that cannot be predicted beforehand: 
a hard process like <i>g g -> g g</i> can develop a shower with a 
<i>g -> b bbar</i> branching, where the <i>b</i> hadronizes to a 
<i>B^0bar</i> that oscillates to a <i>B^0</i> that decays to a 
<i>tau^+</i>. Therefore any change of branching ratios - most of which
are determined from data rather than from first principles anyway -
will not be taken into account in the cross section of a process.
Exceptions exist, but most particles in this class are made to decay
isotropically. Finally, note that all of these particles have a mass 
below 20 GeV.
</li>
</ul>

There is one ambiguous case in this classification, namely the photon.
The <i>gamma^*/Z^0</i> combination contains a low-mass peak when 
produced in a hard process. On the other hand, photons can participate 
in shower evolution, and therefore a photon originally assumed
massless can be assigned an arbitrarily high mass when it is allowed
to branch into a fermion pair. In some cases this could lead to 
doublecounting, e.g. between processes such as 
<i>f fbar -> (gamma^*/Z^0) (gamma^*/Z^0)</i>,
<i>f fbar -> (gamma^*/Z^0) gamma</i> and 
<i>f fbar -> gamma gamma</i>. Here it make sense to limit the
lower mass allowed for the <i>gamma^*/Z^0</i> combination, 
in <code>23:mMin</code> to be the same as the upper limit allowed
for an off-shell photon in the shower evolution, in
<code>TimeShower:mMaxGamma</code>. By default this matching is done 
at 10 GeV.

<p/>
In spite of the above-mentioned differences, the resonances and the 
other particles are all stored in one common 
<?php $filepath = $_GET["filepath"];
echo "<a href='ParticleData.php?filepath=".$filepath."' target='page'>";?>particle data table</a>, so as to offer a 
uniform interface to <?php $filepath = $_GET["filepath"];
echo "<a href='ParticleDataScheme.php?filepath=".$filepath."' target='page'>";?>setting and 
getting</a> properties such as name, mass, charge and decay modes,
also for the <?php $filepath = $_GET["filepath"];
echo "<a href='ParticleProperties.php?filepath=".$filepath."' target='page'>";?>particle properties</a> 
in the event record. Some methods are specific to resonances, however,
in particular for the calculation of partial widths and thereby of
branching ratio. For resonances these can be calculated dynamically,
set up at initialization for the nominal mass and then updated to the
current mass when these are picked according to a Breit-Wigner resonance 
shape. 

<h3>Special properties and methods for resonances</h3>

The method <code>ParticleDataTable::isResonance(id)</code> allows you to 
query whether a given particle species is considered a resonance or not. 
You can also change the default value of this flag in the normal way, 
e.g. <code>pythia.readString("id:isResonance = true")</code>.

<p/>
An option with a a forced width can be set with the 
<code>id:doForceWidth</code> flag as above, and queried with
<code>ParticleDataTable::doForceWidth(id)</code>. It is by default 
<code>off</code>, and should normally 
so remain. If switched <code>on</code> then the width stored in 
<code>id:mWidth</code> is strictly used to describe the Breit-Wigner 
of the resonance. This is unlike the normal behaviour of standard
resonances such as the <i>Z^0</i>, <i>W^+-</i>, <i>t</i> or
<i>h^0</i>, which have explicit decay-widths formulae encoded, 
in classes derived from the <code>ResonanceWidths</code> base class.
These formulae are used, e.g., to derive all the Higgs partial
widths as a function of the Higgs mass you choose, and at initialization
overwrites the existing total width value. The reason for forcing the 
width  to another value specified by you would normally more have to do 
with experimental issues than with physics ones, e.g. how sensitive your
detector would be to changes in the Higgs width by a factor of two.
A warning is that such a rescaling could modify the cross section of
a process correspondingly for some processes, while leaving it 
(essentially) unchanged for others (as would seem most logical), 
depending on how these were encoded.  

<p/>
If a resonance does not have a class of its own, with hardcoded equations 
for all relevant partial widths, then a simpler object will be created
at initialization. This object will take the total width and branching
ratios as is (with the optional variations explained in the next section),
and thus the rescaling approach makes no sense.  

<p/>
Mainly for internal usage, the <code>ParticleDataTable</code> contain
some special methods that are only meaningful for resonances:
<ul>
<li><code>resInit(...)</code> to initialize a resonance, possibly
including a recalculation of the nominal width to match the nominal 
mass;
<li><code>resWidth(...)</code> to calculate the partial and total widths
at the currently selected mass;
<li><code>resWidthChan(...)</code> to return pretabulated widths for
channels (only used for Higgs decays, where some channels require
numerical integration);
<li><code>resOpenFrac(...)</code> to return the fraction of the total 
width that is open by the decay channel selection made by users (based on
the choice of <?php $filepath = $_GET["filepath"];
echo "<a href='ParticleDataScheme.php?filepath=".$filepath."' target='page'>";?><code>onMode</code></a>
for the various decay channels, recursively calculated for sequential 
decays);
<li><code>resWidthRescaleFactor(...)</code> returns the factor by which 
the internally calculated PYTHIA width has to be rescaled to give the
user-enforced width.
</ul>
These methods actually provide an interface to the classes derived from
the <code>ResonanceWidths</code> base class, to describe various 
resonances.  
 
<h3>Modes for Matrix Element Processing</h3>

The <code>meMode()</code> value for a decay mode is used to specify 
<?php $filepath = $_GET["filepath"];
echo "<a href='ParticleDecays.php?filepath=".$filepath."' target='page'>";?>nonisotropic decays or the conversion of 
a parton list into a set of hadrons</a> in some channels of normal 
particles. For resonances it can also take a third function, namely 
to describe how the branching ratios and widths of a resonance should 
be rescaled as a function of the current mass of the decaying resonance. 
The codes are not intended for channels where the 
<code>ResonanceWidths</code> class already contains the correct formulae 
for the various partial widths, but rather for when the only information 
available is the on-shell branching ratios and the total width, 
typically read in from an external program. If no code is set, i.e. 
if one defaults to 0, then no special threshold behaviour is imposed.
If the mother fluctuates down in mass, to below the nominal threshold,
it is assumed that one of the daughters could also fluctuate down to
keep the channel open. (If not, there may be problems later on.)
<ul>
<li>101 : use a step threshold to calculate the current partial width 
from the stored total width and branching ratio, i.e. a channel is 
switched off when the sum of the daughter on-shell masses is above the
current mother mass.</li>
<li>102 : use a smooth threshold factor 
<i>beta = sqrt( (1 - m_1^2/m_2 - m_2^2/m^2)^2 - 4 m_1^2 m_2^2/m^4)</i>
for two-body decays and <i>sqrt(1 - Sum_i m_i / m)</i> for multibody
ones. The current partial width is defined as the product of this 
threshold factor times the stored total width and branching ratio.
Specifically, it is thereby assumed that the stored branching ratio
did not take into account such a factor.</li>
<li>103 : use the same kind of threshold factor as for 102 above,
but assume that such a threshold factor had been used when the default
on-shell total width and branching ratios were calculated, so that one
should additionally divide by the on-shell threshold factor. Specifically,
this will give back the stored branching ratios for on-shell mass,
unlike the 102 option. To avoid division by zero, or in general 
unreasonably big rescaling factors, a lower limit on the value of the
on-shell threshold factor is imposed, hardcoded to 0.1. (In cases
where a big rescaling is intentional, code 102 would anyway be more
appropriate.) </li>
</ul>

</body>
</html>

<!-- Copyright (C) 2007 Torbjorn Sjostrand -->

