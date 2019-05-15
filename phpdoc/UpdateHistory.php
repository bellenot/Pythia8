<html>
<head>
<title>Update History</title>
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

<form method='post' action='UpdateHistory.php'>

<h2>Update History</h2>

These update notes describe major updates relative to the baseline 
PYTHIA 8.100 version. However, they are less extensive than the 
corresponding update notes for PYTHIA 6. There are three main 
reasons for this:

<ul>

<li>The manual contained on these html/php pages is kept up to date.
(However, the "Brief Introduction" may not always be.)</li>

<li>8.1 is a quite new code, so there are many minor changes that, 
if all were to be documented, would hide the key ones.</li>

<li>8.1 is not yet used for "mission critical" applications, 
so there is less need to trace changed behaviour.</li>

</ul>

<h3>Main news by version</h3>

<ul>

<li>8.101: 10 November 2007
<ul>

<li>New option to initialize with arbitrary beam directions<br/>
<code>pythia.init( idA, idB, pxA, pyA, pzA, pxB, pyB, pzB)</code></li>

<li>The <code>LHAevnt</code> and <code>LHAinit</code> classes have been 
joined into a new <code>LHAup</code> one, with new options that allow 
the writing of a Les Houches Event File.</li>

</ul>
</li>

<li>8.102: 6 December 2007
<ul>

<li>Limited capability to use two different <code>Pythia</code> instances 
for signal + pileup event generation, see <code>main19.cc</code> for an 
example.</li> 

<li>Added capability to set <?php $filepath = $_GET["filepath"];
echo "<a href='BeamParameters.php?filepath=".$filepath."' target='page'>";?>beam energy spread 
and beam vertex</a>.
<br/>
<b>Warning:</b> as a consequence, the old <code>Beams</code> group of 
settings has been moved to <code>BeamRemnants</code>, and 
<code>Beams</code> is now instead used for machine beam parameters. 
Therefore also some <code>Main</code> settings of this character have been 
regrouped, as follows:
<table cellspacing="5">
<tr> <td>8.100 setting  </td> <td>has been moved to </td> </tr>
<tr> <td><code>Main:idA          </code></td>
     <td><code>Beams:idA         </code></td> </tr>
<tr> <td><code>Main:idB          </code></td>
     <td><code>Beams:idB         </code></td> </tr>
<tr> <td><code>Main:inCMframe    </code></td>
     <td>(<code>Beams:frameType</code>) </td> </tr>
<tr> <td><code>Main:eCM          </code></td>
     <td><code>Beams:eCM         </code></td> </tr>
<tr> <td><code>Main:eA           </code></td>
     <td><code>Beams:eA          </code></td> </tr>
<tr> <td><code>Main:eB           </code></td>
     <td><code>Beams:eB          </code></td> </tr>
<tr> <td><code>Main:LHEF         </code></td>
     <td><code>Beams:LHEF        </code></td> </tr>
</table></li>

</ul>
</li>

<li>8.103: 22 January 2008
<ul>

<li>Updated HepMC conversion routine.</li>

<li>In the <code>Event</code> class the <code>=</code> and 
<code>=+</code> methods have been overloaded to allow the copying 
or appending of event records. Illustrated in <code>main19.cc</code>.</li>

</ul>
</li>

<li>8.104: 14 February 2008
<ul>

<li>Updated configure scripts.</li>

<li>The <code>SusyLesHouches</code> class updated to handle 
SLHA version 2.</li>

<li>The <code>forceHadronLevel()</code> method introduced for standalone 
hadronization.</li>

<li><code>main15.cc</code> illustrated how either full hadronization or 
only decays of some particles can be looped over for the rest of the 
event retained.</li>

<li>The html and php page formatting improved with 
cascading style sheets.</li>

<li>The static <code>ErrorMsg</code> class has been removed and 
its functionality moved into the non-static <code>Info</code> class,
in the renamed Info file.</li>

</ul>
</li>

<li>8.105: 24 February 2008
<ul>

<li>Further reduction of the use of static, with related code changes.
This should allow to have several almost independent <code>Pythia</code> 
instances. Some static classes still remain, however, notably for
random number generation and particle properties.</li>

<li>Several minor improvements and new options.</li>

</ul>
</li>

<li>8.106: 11 March 2008
<ul>

<li>Improved handling of the Higgs width, relevant for massive and thereby
broad resonance shapes.</li>

</ul>
</li>

<li>8.107: 17 March 2008
<ul>

<li>Correction in the event record, so that the beam particles in line 
1 and 2 do not have any mother according to the <code>motherList</code>
method. Previously the "system" entry in line 0 was counted as their 
mother, which gave rise to an unexpected extra vertex in the conversion 
to the HepMC format.</li>

</ul>
</li>

<li>8.108: 1 May 2008
<ul>

<li>Support for HepMC version 1 is removed, to simplify the code and 
reflect the evolution of the field.</li>

<li>Status codes are stored in HepMC only as 1 for existing and 2 for
decayed or fragmented particles (whereas previously the original PYTHIA
codes were used for the latter).</li>

<li>Parton densities are stored in HepMC as <i>xf(x,Q^2)</i> 
rather than the <i>f(x,Q^2)</i> used in (some) previous versions.</li>

<li>The SusyLesHouches class has ben updated so that reading is fully
compatible with the SLHA2 standard. </li>

<li>The matrix elements for neutralino pair production have now been
completed and checked.</li>

<li>A new compilation option <code>-Wshadow</code> is introduced and 
code is rewritten at all places where this option gave warnings.</li>

<li>Minor library correction to allow compilation with gcc 4.3.0.</li>

<li>Ensure that <i>alpha_strong</i> does not blow up, by introducing 
a minimal scale somewhat above <i>Lambda_3</i> (roughly where
<i>alpha_strong = 10</i>).</li>

<li>New methods <code>isValence1()</code> and <code>isValence2()</code> 
in the <code>Info</code> class.</li>

<li>Protection against division by zero in calculation of decay vertex
(for zero-mass gluons with zero lifetime, where there should be no
displacement).</li>

</ul>
</li>

<li>8.109: 6 May 2008
<ul>

<li>Inclusion of some first code for rescatterng in multiple interactions;
currently only for internal tryout.
</li>

</ul>
</li>

<li>8.110: 18 June 2008
<ul>

<li>Freeze GRV 94L distribution at small <i>Q^2</i> to avoid blowup.</li>

<li>Further evolution of rescattering to include colour reconnection 
and more, but not yet FSR or ISR; still only for internal tryout.</li>

</ul>
</li>

<li>8.111: 9 July 2008
<ul>

<li>Further evolution of rescattering to include FSR (with some issues 
still unsolved) but not yet ISR; still only for internal tryout.
</li>

</ul>
</li>

<li>8.112: 12 September 2008
<ul>

<li>The <code>pythia.readFile(...)</code> method can now alternatively take 
an <code>istream</code> as argument instead of a <code>filename</code>.</li>

<li>Minor bug correction in <code>PartonLevel.cc</code>; the bug could 
(rarely) give a segmentation fault.</li>

<li>Resolution of several (but not all) rescattering issues for FSR; still 
only for internal tryout.</li>

</ul>
</li>

<li>8.113: 2 October 2008
<ul>

<li>A new class <code>PartonSystems</code> is introduced to 
keep track of which partons in the event record belong to which 
subcollision system, plus some further information on each subsystem.
It takes over functionality previously found as part of the 
<code>Event</code> class, but leaves room for future growth.</li>

<li>Development of rescattering in combination with ISR; still only for 
internal tryout.</li>

</ul>
</li>

<li>8.114: 22 October 2008
<ul>

<li>Add capability to link to FastJet, with <code>main61.cc</code> 
as new example.</li>

<li>Add optional model, wherein an increased <i>pT0</i> turnoff scale 
for MI and ISR is used for above-average active events, i.e. events that 
already have several MI's or ISR emissions.</li>

<li>Rescattering description operational (but still experimental) for case 
that FSR is not interleaved, but saved until after MI, ISR and beam remnants 
have been handled.</li>

</ul>
</li>

</ul>

</body>
</html>

<!-- Copyright (C) 2008 Torbjorn Sjostrand -->
