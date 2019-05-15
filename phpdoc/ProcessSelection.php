<html>
<head>
<title>Process Selection</title>
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

<form method='post' action='ProcessSelection.php'>

<h2>Process Selection</h2>

By default all processes are switched off. You should switch on 
those you want to simulate. This may be done at two levels, either
for each individual process or for a group of processes. That is,
a process is going to be generated either if its own flag or its 
group flag is on. There is no built-in construction to switch on
a group and then switch off a few of its members.

<p/>
Each process is assigned an integer code. This code is not used in
the internal administration of events (so having the same code for
two completely different processes would not be a problem), but only 
intended to allow a simpler user separation of different processes. 
Also the process name is available, as a string.

<p/>
To ease navigation, the list of processes has been split into several
separate pages, by main topic. The classification is hopefully
intuitive, but by no means unambiguous. For instance, essentially 
all processes involve QCD, so the "QCD processes" are the ones that
only involve QCD. (And also that is not completely true, once one 
includes all that may happen in multiple interactions.) On these 
separate pages also appear the settings that are completely local
to that particular process class, but not the ones that have a 
broader usage.

<h3><?php $filepath = $_GET["filepath"];
echo "<a href='QCDProcesses.php?filepath=".$filepath."' target='page'>";?>QCD Processes</a></h3>

QCD processes fall in two main categories: soft and hard. The soft ones
contain elastic, diffractive and "minimum-bias" events, together
covering the total cross section. Hard processea are the normal
<i>2 -> 2</i> ones, including charm and bottom production.  

<h3><?php $filepath = $_GET["filepath"];
echo "<a href='ElectroweakProcesses.php?filepath=".$filepath."' target='page'>";?>Electroweak Processes</a></h3>

Prompt-photon, <i>gamma^*/Z^0</i> and <i>W^+-</i> production, 
plus a few processes with <i>t</i>-channel boson exchange.

<h3><?php $filepath = $_GET["filepath"];
echo "<a href='OniaProcesses.php?filepath=".$filepath."' target='page'>";?>Onia Processes</a></h3>

Colour singlet and octet production of charmonium and bottomonium.

<h3><?php $filepath = $_GET["filepath"];
echo "<a href='TopProcesses.php?filepath=".$filepath."' target='page'>";?>Top Processes</a></h3>

Top production, singly or doubly.

<h3><?php $filepath = $_GET["filepath"];
echo "<a href='HiggsProcesses.php?filepath=".$filepath."' target='page'>";?>Higgs Processes</a></h3>

Higgs production, currently within the Standard Model but eventually
also beyond it.

<h3><?php $filepath = $_GET["filepath"];
echo "<a href='SUSYProcesses.php?filepath=".$filepath."' target='page'>";?>SUSY Processes</a></h3>

Production of supersymmetric particles, currently barely begun

<h3>Resonance Decays and Cross Sections</h3>

In addition to the switches and parameters in the process machinery
there also exists the possibility to set the allowed decay channels
of resonances, as explained in the 
<?php $filepath = $_GET["filepath"];
echo "<a href='ParticleDataScheme.php?filepath=".$filepath."' target='page'>";?>Particle Data Scheme</a> description.
For instance, if you study the process <i>q qbar -> H^0 Z^0</i>
you could specify that the <i>Z^0</i> should decay only to
lepton pairs, the <i>H^0</i> only to <i>W^+ W^-</i>, the 
<i>W^+</i> only to a muon and a neutrino, while the <i>W^-</i>
can decay to anything. Unfortunately there are limits to the 
flexibility: you cannot set a resonance to have different properties
in different places of a process, e.g. if instead 
<i>H^0 -> Z^0 Z^0</i> in the above process then the three 
<i>Z^0</i>'s would all obey the same rules.

<p/>
The restrictions on the allowed final states of a process is directly
reflected in the cross section of it. That is, if some final states
are excluded then the cross section is reduced accordingly. Such 
restrictions are built up recursively in cases of sequential decay 
chains. The restrictions are also reflected in the compositions of
those events that actually do get to be generated. For instance,
the relative rates of <i>H^0 -> W^+ W^-</i> and 
<i>H^0 -> Z^0 Z^0</i> are shifted when the allowed sets of 
<i>W^+-</i> and <i>Z^0</i> decay channels are changed.

<p/>
There is one important restriction, however: only those particles that
Pythia treat as resonances enjoy this property. This includes the
<i>W^+-</i>, <i>gamma^*/Z^0</i>, <i>t/tbar</i>, <i>H^0</i> and 
other heavy unstable particles in scenarios of Beyond-the-Standard-Model 
physics. These particles are generated and let to decay as part of the 
"process level" processing, which is where cross sections are handled.
It does <i>not</i> concern particle that are produced and/or decay at
later stages, such as <i>B</i> mesons or <i>tau</i> leptons, or
photons that branch as part of the shower evolution. There simply
would be no way consistently to include the proper bias that should go 
with changed branching ratios. For instance, if you only are interested
in a specific <i>tau</i> decay channel, this <i>tau</i> could come 
from the decay of a <i>B</i> meson that came from a <i>b</i> quark 
produced in the shower evolution, <i>g -> b bbar</i>, and thus be 
many steps removed from the hard process itself.  

</body>
</html>

<!-- Copyright (C) 2007 Torbjorn Sjostrand -->

