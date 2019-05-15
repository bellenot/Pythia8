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
echo "<a href='TopProcesses.php?filepath=".$filepath."' target='page'>";?>Top Processes</a></h3>

Top production, singly or doubly.

<h3><?php $filepath = $_GET["filepath"];
echo "<a href='OniaProcesses.php?filepath=".$filepath."' target='page'>";?>Onia Processes</a></h3>

Colour singlet and octet production of charmonium and bottomonium.

<h3><?php $filepath = $_GET["filepath"];
echo "<a href='SUSYProcesses.php?filepath=".$filepath."' target='page'>";?>SUSY Processes</a></h3>

Production of supersymmetric particles.

</body>
</html>

<!-- Copyright C 2007 Torbjorn Sjostrand -->

