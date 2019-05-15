<html>
<head>
<title>Event Information</title>
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

<form method='post' action='EventInformation.php'>

<h2>Event Information</h2>

The <code>Info</code> class collects various one-of-a-kind information, 
some relevant for all events and others for the current event. 
An object <code>info</code> is a public member of the <code>Pythia</code>
class, so if you e.g. have declared <code>Pythia pythia</code>, the
<code>Info</code> methods can be accessed by 
<code>pythia.info.method()</code>. Most of this is information that 
could also be obtained e.g. from the event record, but is here more
directly available. It is primarily intended for processes generated 
internally in PYTHIA, but many of the methods would work also for
events fed in via the Les Houches Accord.

<h3>List information</h3>

<p/><strong>void Info::list() &nbsp;</strong> <br/>
a listing of most of the information set for the current event. 
  

<h3>The beams</h3>

<p/><strong>int Info::idA() &nbsp;</strong> <br/>
  
<strong>int Info::idB() &nbsp;</strong> <br/>
the identities of the two beam particles. 
  

<p/><strong>double Info::pzA() &nbsp;</strong> <br/>
  
<strong>double Info::pzB() &nbsp;</strong> <br/>
the longitudinal momenta of the two beam particles.
  

<p/><strong>double Info::eA() &nbsp;</strong> <br/>
  
<strong>double Info::eB() &nbsp;</strong> <br/>
the energies of the two beam particles.
  

<p/><strong>double Info::mA() &nbsp;</strong> <br/>
  
<strong>double Info::mB() &nbsp;</strong> <br/>
the masses of the two beam particles.
  

<p/><strong>double Info::eCM() &nbsp;</strong> <br/>
  
<strong>double Info::s() &nbsp;</strong> <br/>
the cm energy and its square for the two beams. 
  

<h3>The event type</h3>

<p/><strong>string Info::name() &nbsp;</strong> <br/>
  
<strong>int Info::code() &nbsp;</strong> <br/>
the name and code of the process that occured.
  

<p/><strong>int Info::nFinal() &nbsp;</strong> <br/>
the number of final-state partons in the hard process.
  

<p/><strong>bool Info::isResolved() &nbsp;</strong> <br/>
are beam particles resolved, i.e. were PDF's used for the process?
  

<p/><strong>bool Info::isDiffractiveA() &nbsp;</strong> <br/>
  
<strong>bool Info::isDiffractiveB() &nbsp;</strong> <br/>
is either beam diffractively excited?
  

<p/><strong>bool Info::isMinBias() &nbsp;</strong> <br/>
is the process a minimum-bias one?
  

<p/><strong>bool Info::isLHA() &nbsp;</strong> <br/>
has the process been generated from external Les Houches Accord 
information?
  

<p/><strong>bool Info::atEndOfFile() &nbsp;</strong> <br/>
true if a linked Les Houches class refuses to return any further 
events, presumably because it has reached the end of the file from 
which events have been read in.
  

<p/><strong>bool Info::hasSub() &nbsp;</strong> <br/>
does the process have a subprocess classification?
Currently only true for minbias and Les Houches events, where it allows 
the hardest collision to be identified. 
  

<p/><strong>string Info::nameSub() &nbsp;</strong> <br/>
  
<strong>int Info::codeSub() &nbsp;</strong> <br/>
  
<strong>int Info::nFinalSub() &nbsp;</strong> <br/>
the name, code and number of final-state partons in the subprocess
that occured when <code>hasSub()</code> is true. For a minimum-bias event 
the <code>code</code> would always be 101, while <code>codeSub()</code> 
would vary depending on the actual hardest interaction, e.g. 111 for 
<i>g g -> g g</i>. For a Les Houches event the <code>code</code> would 
always be 9999, while <code>codeSub()</code> would be the external 
user-defined classification code. The methods below would also provide 
information for such particular subcollisions.  
  

<h3>Hard process parton densities and scales</h3>

<p/><strong>int Info::id1() &nbsp;</strong> <br/>
  
<strong>int Info::id2() &nbsp;</strong> <br/>
the identities of the two partons coming in to the hard process.
  

<p/><strong>double Info::x1() &nbsp;</strong> <br/>
  
<strong>double Info::x2() &nbsp;</strong> <br/>
<i>x</i> fractions of the two partons coming in to the hard process.
  

<p/><strong>double Info::y() &nbsp;</strong> <br/>
  
<strong>double Info::tau() &nbsp;</strong> <br/>
rapidity and scaled mass-squared of the hard-process subsystem, as 
defined by the above <i>x</i> values. 
  

<p/><strong>double Info::pdf1() &nbsp;</strong> <br/>
  
<strong>double Info::pdf2() &nbsp;</strong> <br/>
parton densities <i>x*f(x,Q^2</i> )evaluated for the two incoming 
partons; could be used e.g. for reweighting purposes. 
  

<p/><strong>double Info::QFac() &nbsp;</strong> <br/>
  
<strong>double Info::Q2Fac() &nbsp;</strong> <br/>
the <i>Q^2</i> or <i>Q^2</i> factorization scale at which the 
densities were evaluated.
  

<p/><strong>bool Info::isValence1() &nbsp;</strong> <br/>
  
<strong>bool Info::isValence2() &nbsp;</strong> <br/>
<code>true</code> if the two hard incoming partons have been picked 
to belong to the valence piece of the parton-density distribution, 
else <code>false</code>. Should be interpreted with caution.
Information is not set if you switch off parton-level processing. 
  

<p/><strong>double Info::alphaS() &nbsp;</strong> <br/>
  
<strong>double Info::alphaEM() &nbsp;</strong> <br/>
the <i>alpha_strong</i> and <i>alpha_electromagnetic</i> values used 
for the hard process.
  

<p/><strong>double Info::QRen() &nbsp;</strong> <br/>
  
<strong>double Info::Q2Ren() &nbsp;</strong> <br/>
the <i>Q</i> or <i>Q^2</i> renormalization scale at which 
<i>alpha_strong</i> and <i>alpha_electromagnetic</i> were evaluated.
  

<h3>Hard process kinematics</h3>

<p/><strong>double Info::mHat() &nbsp;</strong> <br/>
  
<strong>double Info::sHat() &nbsp;</strong> <br/>
the invariant mass and its square for the hard process.
  

<p/><strong>double Info::tHat() &nbsp;</strong> <br/>
  
<strong>double Info::uHat() &nbsp;</strong> <br/>
the remaining two Mandelstam variables; only defined for <i>2 -> 2</i>
processes. 
  

<p/><strong>double Info::pTHat() &nbsp;</strong> <br/>
  
<strong>double Info::pT2Hat() &nbsp;</strong> <br/>
transverse momentum and its square in the rest frame of a <i>2 -> 2</i>
processes. 
  

<p/><strong>double Info::m3Hat() &nbsp;</strong> <br/>
  
<strong>double Info::m4Hat() &nbsp;</strong> <br/>
the masses of the two outgoing particles in a <i>2 -> 2</i> processes. 
  

<p/><strong>double Info::thetaHat() &nbsp;</strong> <br/>
  
<strong>double Info::phiHat() &nbsp;</strong> <br/>
the polar and azimuthal scattering angles in the rest frame of 
a <i>2 -> 2</i> process.
  

<h3>Event weight and activity</h3>

<p/><strong>double Info::weight() &nbsp;</strong> <br/>
weight assigned to the current event. Is normally 1 and thus uninteresting. 
However, for Les Houches events some strategies allow negative weights, 
which then after unweighting lead to events with weight -1. There are also 
strategies where no unweighting is done, and therefore a nontrivial event 
weight must be used e.g. when filling histograms. 
  

<p/><strong>int Info::nISR() &nbsp;</strong> <br/>
  
<strong>int Info::nFSRinProc() &nbsp;</strong> <br/>
  
<strong>int Info::nFSRinRes() &nbsp;</strong> <br/>
the number of emissions in the initial-state showering, in the final-state
showering excluding resonance decys, and in the final-state showering
inside resonance decays, respectively.
  

<p/><strong>double Info::pTmaxMI() &nbsp;</strong> <br/>
  
<strong>double Info::pTmaxISR() &nbsp;</strong> <br/>
  
<strong>double Info::pTmaxFSR() &nbsp;</strong> <br/>
Maximum <i>pT</i> scales set for MI, ISR and FSR, given the 
process type and scale choice for the hard interactions. The actual
evolution will run down from these scales.
  

<h3>Multiple interactions</h3>

<p/><strong>double Info::bMI() &nbsp;</strong> <br/>
the impact parameter <i>b</i> assumed for the current collision when
multiple interactions are simulated. Is not expressed in any physical
size (like fm), but only rescaled so that the average should be unity 
for minimum-bias events (meaning less than that for events with hard
processes). 
  

<p/><strong>double Info::enhanceMI() &nbsp;</strong> <br/>
The choice of impact parameter implies an enhancement or depletion of
the rate of subsequent interactions, as given by this number. Again
the average is normalized be unity for minimum-bias events (meaning 
more than that for events with hard processes).  
  

<p/><strong>int Info::nMI() &nbsp;</strong> <br/>
the number of hard interactions in the current event. Is 0 for elastic
and diffractive events, and else at least 1, with more possible from
multiple interactions.
  

<p/><strong>int Info::codeMI(int i) &nbsp;</strong> <br/>
  
<strong>double Info::pTMI(int i) &nbsp;</strong> <br/>
the process code and transverse momentum of the <code>i</code>'th 
subprocess, with <code>i</code> in the range from 0 to
<code>nMI() - 1</code>. The values for subprocess 0 is redundant with
information already provided above.  
  

<p/><strong>int Info::iAMI(i) &nbsp;</strong> <br/>
  
<strong>int Info::iBMI(i) &nbsp;</strong> <br/>
are normally zero. However, if the <code>i</code>'th subprocess is
a rescattering, i.e. either or both incoming partons come from the 
outgoing state of previous scatterings, they give the position in the
event record of the outgoing-state parton that rescatters. 
<code>iAMI</code> and <code>iBMI</code> then denote partons coming from 
the first or second beam, respectively.
  

<h3>Cross sections</h3>

Here are the currently available methods related to the event sample 
as a whole. While continuously updated during the run, it is recommended
only to study these properties at the end of the event generation, 
when the full statistics is available.

<p/><strong>long Info::nTried() &nbsp;</strong> <br/>
  
<strong>long Info::nSelected() &nbsp;</strong> <br/>
  
<strong>long Info::nAccepted() &nbsp;</strong> <br/>
the total number of tried phase-space points, selected hard processes
and finally accepted events, summed over all allowed subprocesses.
The first number is only intended for a study of the phase-space selection
efficiency. The last two numbers usually only disagree if the user introduces 
some veto during the event-generation process; then the former is the number 
of acceptable events found by PYTHIA and the latter the number that also
were approved by the user. If you set <?php $filepath = $_GET["filepath"];
echo "<a href='ASecondHardProcess.php?filepath=".$filepath."' target='page'>";?>a 
second hard process</a> there may also be a mismatch. 
  

<p/><strong>double Info::sigmaGen() &nbsp;</strong> <br/>
  
<strong>double Info::sigmaErr() &nbsp;</strong> <br/>
the estimated cross section and its estimated error,
summed over all allowed subprocesses, in units of mb. The numbers refer to
the accepted event sample above, i.e. after any user veto. 
  

</body>
</html>

<!-- Copyright (C) 2009 Torbjorn Sjostrand -->
