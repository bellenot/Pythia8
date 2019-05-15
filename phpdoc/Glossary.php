<html>
<head>
<title>Glossary</title>
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

<form method='post' action='Glossary.php'>

<h2>Glossary</h2>

<dl>

<dt>BR</dt>
<dd>Beam Remnants; not much used since it may be confused with 
Branching Ratio</dd>

<dt>FSR</dt>
<dd>Final-State Radiation, implemented in terms of timelike showers</dd>

<dt>ISR</dt>
<dd>Initial-State Radiation, implemented in terms of spacelike showers</dd>

<dt>MI</dt>
<dd>Multiple Interactions, i.e. several (more or less) independent 
parton-parton subcollisions as part of a hadron-hadron event (sometimes
also called MPI, with P for parton or parton-parton)</dd>

<dt>PDF</dt>
<dd>Parton Distribution Function (alternatively Parton Density 
Function)</dd>

<dt>pileup</dt>
<dd>several hadron-hadron collisions in a bunch crossing</dd>

<dt>RPP</dt>
<dd>Review of Particle Physics, the biannual review by the ParticleData Group
(PDG) from which many Standard-Model parameter values and much particle data 
has been taken (but, given the poor data on many hadron resonances, a lot of
extra (guess)work is needed)</dd>

<dt>setting</dt>
<dd>collectively used to denote all the boolean <code>flag</code>, 
integer <code>mode</code>, double-precision <code>parm</code>
and string <code>word</code> variables that can be set by the user 
to steer the behaviour of a run; normally particle data are considered
separately but clearly are closely related</dd>

</dl>

</body>
</html>

<!-- Copyright C 2007 Torbjorn Sjostrand -->
