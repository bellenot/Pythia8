<html>
<head>
<title>Program Files</title>
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

<form method='post' action='ProgramFiles.php'>

<h2>Program Files</h2>

The code is subdivided into a set of files, mainly by physics 
task. Each file typically contains one main class, but often
with a few related helper classes that are not used elsewhere in 
the program. Normally the files come in pairs.
<ul>
<li>A header file, <code>.h</code> in the <code>include</code> 
subdirectory, where the public interface of the class is declared, 
and inline methods are defined.</li>
<li>A source code file, <code>.cc</code> in the <code>src</code> 
subdirectory, where most of the methods are implemented.</li>
</ul>
During compilation, related dependency files, <code>.d</code>, and 
compiled code, <code>.o</code> are created in the <code>tmp</code>
subdirectory.

<p/>
In part the <code>.xml</code> documentation files in the 
<code>xmldoc</code> subdirectory have matching names, but the match 
is broken by the desire to group topics more by user interaction than 
internal operation. On these pages the function of the different code 
files is summarized. Currently, each <code>.xml</code> file is also 
translated into an <code>.html</code> one in the 
<code>htmldoc</code> subdirectory, to allow easy viewing of the 
contents in a web browser, and an <code>.php</code> one in 
<code>phpdoc</code>, for more sophisticated interactivity
if that subdirectory is installed on a web server.

<p/><code>file&nbsp; </code><strong> Pythia &nbsp;</strong> <br/>
is the main class, that administrates the whole event generation 
process by making use of all the others classes. Objects of most
other classes reside (directly or indirectly) inside <code>Pythia</code>, 
so only a <code>Pythia</code> object needs to be explicitly instantiated 
and addressed   by the user.
  

<p/><code>file&nbsp; </code><strong> Event &nbsp;</strong> <br/>
contains the event record, which basically is a vector of particles.
This file also contains the <code>Particle</code> class, used by 
<code>Event</code>. <code>Pythia</code> uses two <code>Event</code> 
objects, one for the process-level record (<code>process</code>) and 
one for the complete event (<code>event</code>).
  

<p/><code>file&nbsp; </code><strong> Information &nbsp;</strong> <br/>
is a simple container that gives access to some information on the
nature of the current process, such as Mandelstam variables.
Also contains a small database for errors and warnings encountered 
during program execution.
  

<p/><code>file&nbsp; </code><strong> Analysis &nbsp;</strong> <br/>
contains routines to analyze events. Currently it can do sphericity, 
thrust, Lund/Jade/Durham jet clustering, and cone-jet finding.
  

<p/><code>file&nbsp; </code><strong> ProcessLevel &nbsp;</strong> <br/>
handles the generation of the (hard) process that sets the character 
of the event. This involves either using internally implemented 
processes  (few so far), using a runtime interface to the Fortran 
PYTHIA 6 process library, or using routines for reading in 
LHA events. 
  

<p/><code>file&nbsp; </code><strong> PartonLevel &nbsp;</strong> <br/>
turns the (hard) process above into a complete set of partons, by 
adding initial- and final-state radiation, multiple parton--parton
interactions, and beam remnants.
  

<p/><code>file&nbsp; </code><strong> HadronLevel &nbsp;</strong> <br/>
turns the parton-level event above into a set of outgoing particles,
by applying string fragmentation (with special treatment for low-mass
systems) and secondary decays.
  

<p/><code>file&nbsp; </code><strong> ProcessContainer &nbsp;</strong> <br/>
packages the information on a given subprocess, combining the 
phase-space selection and cross-section evaluation machineries
with some statistics information. Also sets up the list of processes
to be studied in a run.  
  

<p/><code>file&nbsp; </code><strong> SigmaTotal &nbsp;</strong> <br/>
contains parametrizations of total, elastic and diffractive hadronic 
cross sections.
  

<p/><code>file&nbsp; </code><strong> SigmaProcess &nbsp;</strong> <br/>
contains the base class and derived classes for the evaluation of
different matrix elements. In order to keep this file from becoming 
too big, actual cross sections are found in the following four files
of derived classes.
  

<p/><code>file&nbsp; </code><strong> SigmaQCD &nbsp;</strong> <br/>
contains the cross sections and matrix elements for soft and hard
QCD processes.
  

<p/><code>file&nbsp; </code><strong> SigmaEW &nbsp;</strong> <br/>
contains the cross sections and matrix elements for electroweak
processes involving photons, <i>Z^0</i>'s and <i>W^+-</i>'s .
  

<p/><code>file&nbsp; </code><strong> SigmaOnia &nbsp;</strong> <br/>
contains the cross sections and matrix elements charmonium and 
bottomonium production.
  

<p/><code>file&nbsp; </code><strong> SigmaSUSY &nbsp;</strong> <br/>
contains the cross sections and matrix elements for Supersymmetric
processes.
  

<p/><code>file&nbsp; </code><strong> SusyLesHouches &nbsp;</strong> <br/>
contains information on SUSY parameters and particle data as specified
by the SUSY Les Houches Accord.
  

<p/><code>file&nbsp; </code><strong> InFlux &nbsp;</strong> <br/>
keeps track of the allowed combinations of incoming flavours for a
process, and stores the in-state-specific factors in a cross section.
This in particular means the parton densities, but also charges,
CKM matrix elements, and so on.
  

<p/><code>file&nbsp; </code><strong> PhaseSpace &nbsp;</strong> <br/>
selects a point in phase space for the hard-process generation,
optimized to give improved Monte Carlo efficiency.
  

<p/><code>file&nbsp; </code><strong> PartonDistributions &nbsp;</strong> <br/>
contains parton distribution functions for the proton and electron. 
Currently very simple, with only two <i>p</i> parametrizations 
and one <i>e</i> ditto available, but it is possible to link in 
external sets.
  

<p/><code>file&nbsp; </code><strong> LHAPDFInterface &nbsp;</strong> <br/>
is a header file only, with interfaces to the key LHAPDF routines,
as needed for a runtime interface. There is a file
<code>lhapdfdummy/LHAPDFdummy.cc</code> with matching dummy
implementations, however. This file is used to build a separate 
<code>liblhapdfdummy</code> library, to be linked when the LHAPDF 
library is not used, so as to avoid problems with undefined references. 
  

<p/><code>file&nbsp; </code><strong> BeamParticle &nbsp;</strong> <br/>
contains information on all partons extracted from one of the two
beams. Defines modified parton distributions accordingly during the
showering and multiple interactions processing, thereby extending on 
the one-particle-inclusive distributions defined by the previous class.
Finds the internal structure for a beam remnant.
  

<p/><code>file&nbsp; </code><strong> LesHouches &nbsp;</strong> <br/>
gives the possibility to feed in parton configurations for the 
subsequent event generation. Two containers are defined, one for 
initialization and one for events, that can be read from 
<code>Pythia</code>. Should be linked to external programs or files.
  

<p/><code>file&nbsp; </code><strong> LHAFortran &nbsp;</strong> <br/>
is a header file only, for two classes derived from the above LesHouches 
ones, to be used for runtime interfacing to Fortran programs, such as
PYTHIA 6. 
  

<p/><code>file&nbsp; </code><strong> Pythia6Interface &nbsp;</strong> <br/>
is a header file only, with interfaces to the key PYTHIA 6 routines,
as needed for a runtime interface. 
  

<p/><code>file&nbsp; </code><strong> ResonanceDecays &nbsp;</strong> <br/>
decays resonances as part of the hard-process stage, eventually including
angular correlations between decay products.
  

<p/><code>file&nbsp; </code><strong> ResonanceProperties &nbsp;</strong> <br/>
encodes some properties of resonances, such as the dynamic calculation 
of widths.
  

<p/><code>file&nbsp; </code><strong> TimeShower &nbsp;</strong> <br/>
performs timelike final-state transverse-momentum-ordered 
shower evolution.
  

<p/><code>file&nbsp; </code><strong> SpaceShower &nbsp;</strong> <br/>
performs spacelike initial-state transverse-momentum-ordered 
shower evolution.
  

<p/><code>file&nbsp; </code><strong> MultipleInteractions &nbsp;</strong> <br/>
performs multiple parton-parton interactions.
  

<p/><code>file&nbsp; </code><strong> BeamRemnants &nbsp;</strong> <br/>
adds primordial <i>kT</i> to the set of hard subsystems,
and combines these subsystems with the two beam remnants to provide 
the overall energy-momentum picture. Also ties together all the 
colour lines into consistent singlets. 
  

<p/><code>file&nbsp; </code><strong> StringFragmentation &nbsp;</strong> <br/>
performs string fragmentation of a given set of partons.
  

<p/><code>file&nbsp; </code><strong> MiniStringFragmentation &nbsp;</strong> <br/>
performs string fragmentation in cases where the colour singlet
subsystem mass is so small that one or at most two primary hadrons 
should be produced from it.
  

<p/><code>file&nbsp; </code><strong> FragmentationFlavZpT &nbsp;</strong> <br/>
contains the classes for describing the fragmentation steps in
flavour and in longitudinal and transverse momentum.
  

<p/><code>file&nbsp; </code><strong> FragmentationSystems &nbsp;</strong> <br/>
defines some containers of parton systems, for use in 
the fragmentation routines.
  

<p/><code>file&nbsp; </code><strong> ParticleDecays &nbsp;</strong> <br/>
performs the decays of all normal unstable hadrons and leptons, i.e. 
in mass up to and including <i>b bbar</i> systems. It is not 
intended for decays of electroweak resonances, like <i>Z^0</i>.  
  

<p/><code>file&nbsp; </code><strong> StandardModel &nbsp;</strong> <br/>
contains the running <i>alpha_strong</i>, with <i>Lambda</i> 
matching at flavour thresholds, the running <i>alpha_em</i>,
CKM mixing matrices, and a few other parameters such as
<i>sin^2(theta_W)</i>.
  

<p/><code>file&nbsp; </code><strong> UserHooks &nbsp;</strong> <br/>
Provides a way for a user to study the event at a few intermediate
stages of evolution, to reject the event as a whole or to modify
it cross-section weight. 
  

<p/><code>file&nbsp; </code><strong> HepMCInterface &nbsp;</strong> <br/>
contains an interface to convert the PYTHIA 8 event record into the 
HepMC format. The <code>HepMCInterface.cc</code> file is located in
the subdirectory <code>hepmcinterface</code> and is used to build a 
separate <code>libhepmcinterface</code> library. 
  

<p/><code>file&nbsp; </code><strong> Settings &nbsp;</strong> <br/>
contains a database of all flags, modes and parameters that determine 
the performance of the generator. Initial values are set from the 
contents of the <code>.xml</code> files, but these values can then 
be changed by the user.
  

<p/><code>file&nbsp; </code><strong> ParticleData &nbsp;</strong> <br/>
contains a database of all necessary particle data (masses, names, ..)
and decay channels.
  

<p/><code>file&nbsp; </code><strong> Basics &nbsp;</strong> <br/>
contains some basic facilities of general use: a random number
generator <code>Rndm</code>, a four-vector class <code>Vec4</code>, and a 
histogram class <code>Hist</code>.  
  

<p/><code>file&nbsp; </code><strong> PythiaStdlib &nbsp;</strong> <br/>
is only a <code>.h</code> file, containing most of the <code>Stdlib</code> 
headers used in Pythia 8, with <code>using</code> directives. Also 
a few simple inline methods.
  

<p/><code>file&nbsp; </code><strong> PythiaComplex &nbsp;</strong> <br/>
is only a <code>.h</code> file, containing a <code>typedef</code> for 
double precision complex numbers.
  

</body>
</html>

<!-- Copyright C 2007 Torbjorn Sjostrand -->
