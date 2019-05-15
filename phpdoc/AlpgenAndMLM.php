
<html>
<head>
<title>ALPGEN and MLM merging</title>
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

<form method='post' action='AlpgenAndMLM.php'>

<h2>ALPGEN and MLM merging</h2>

This manual page describes the ALPGEN [<a href="Bibliography.php" target="page">Man03</a>] and MLM merging
[<a href="Bibliography.php" target="page">Man02</a>] interfaces for PYTHIA 8. While future versions of
ALPGEN will be able to write out events in LHEF format, previous
versions always output events in an ALPGEN native format (a combination
of a ".unw" and a "_unw.par" file). The ALPGEN component of this code
contains a reader for this native format (for unweighted events), as
well as parameter reading for both ALPGEN native and LHE file formats.
Although designed to work together with the MLM component, they are
implemented entirely separately and it is possible to use one without
the other.

<p/>
It should be noted that all the functionality described here is provided
through external routines, and therefore the presence of these features is
dependent on the main program being used. This structure allows for the
easy extensibility of the merging scheme. The files of interest are located
in the <code>examples/</code> subdirectory:
<ul>
<li>
<code>MLMhooks.h</code> : provides a
<?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?><code>UserHooks</code></a> derived
class for performing the MLM merging procedure. All <code>MLM:*</code>
options, described below, are implemented in this file. Further
technical details of this class are given at the end of this manual
page.
</li>
<li>
<code>LHAupAlpgen.h</code> : provides three classes for the reading of
ALPGEN event and parameter files. <code>LHAupAlpgen</code> is an 
<?php $filepath = $_GET["filepath"];
echo "<a href='LesHouchesAccord.php?filepath=".$filepath."' target='page'>";?><code>LHAup</code></a> derived
class for reading in ALPGEN native format event files.
<code>AlpgenPar</code> is a class for the parsing of ALPGEN parameter
files, making the information available through a simple interface.
<code>AlpgenHooks</code> is a
<?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?><code>UserHooks</code></a> derived class that
provides the <code>Alpgen:*</code> options, described below. Further
technical details of these classes are given at the end of this manual
page.
</li>
<li>
<code>main32.cc, main32.cmnd</code> : a sample main program and card
file showing the usage of the previous two files. It combines the Alpgen
and MLM UserHooks classes together, such that the functionality of both
is available, and reads in a sample ALPGEN event file while performing
the MLM merging procedure. Some commented out sets of options are
provided in the card file, which can be activated to try different
merging setups.
</li>
<li>
<code>main32.unw, main32_unw.par</code> : an ALPGEN format event and
parameter file containing 100 W + 3 jet events. It is not feasible
to package large event files with the PYTHIA distribution, but this
sample is enough to show the different components in action.
</li>
</ul>

<h2>ALPGEN</h2>

<b>NB: these options are provided by the AlpgenHooks class,
which must be loaded for this functionality to be present</b>

<p/>
ALPGEN event files that have been written out in LHEF format should be
read in through the normal LHEF machinery
(see <?php $filepath = $_GET["filepath"];
echo "<a href='BeamParameters.php?filepath=".$filepath."' target='page'>";?>beam parameters</a>). Files in
ALPGEN's native format, instead, may be processed using the
<code>Alpgen:file</code> option below. When using this option, the
ALPGEN parameter file is stored in the PYTHIA Info object under the key
<code>AlpgenPar</code>, see the "Header information" section of the
<?php $filepath = $_GET["filepath"];
echo "<a href='EventInformation.php?filepath=".$filepath."' target='page'>";?>Event Information</a> manual page for
more details. Processes not implemented by the PYTHIA 6 interface
supplied with ALPGEN are also not implemented here.

<p/>
When reading in ALPGEN native event files, some momenta are shifted by
the file reader to ensure energy-momentum conservation. The magnitude of
these shifts should be small (around the MeV level in the worst case)
and warnings will be produced if they are above a set threshold. A large
number of warnings may signify unexpected behaviour and should
potentially be investigated. It is also known that certain event
classes, for example an event with both light and heavy <i>b</i>
quarks may give rise to these warnings.

<p/>
The ALPGEN file reader supports the reading of the event and parameter
files in gzip format with file extensions ".unw.gz" and "_unw.par.gz"
respectively. This requires the use of external libraries, however, and
the <code>README</code> file in the main directory contains instructions
on how to enable this.

<p/>
All other <code>Alpgen:*</code> options apply to both LHE and native
file formats, and include options to guide the MLM merging procedure
based on the parameters that are read in with the events file.

<br/><br/><table><tr><td><strong>Alpgen:file  </td><td></td><td> <input type="text" name="1" value="void" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>void</strong></code>)</td></tr></table>
This option is used to read in ALPGEN format event files. Using this option
overrides any previously set beam options inside PYTHIA. The path to the
files, not including any file extension, should be provided e.g. for input
files <i>input_unw.par</i> and <i>input.unw</i>, the value
<i>input</i> should be used.
  

<br/><br/><strong>Alpgen:setMasses</strong>  <input type="radio" name="2" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="2" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
When switched on, any particle masses provided by ALPGEN are set in
the PYTHIA <?php $filepath = $_GET["filepath"];
echo "<a href='ParticleDataScheme.php?filepath=".$filepath."' target='page'>";?>particle database</a>.
  

<br/><br/><strong>Alpgen:setMLM</strong>  <input type="radio" name="3" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="3" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
When switched on, the merging parameters (see below) are set according to
the ALPGEN hard process cuts:
<ul>
<li> <i>MLM:eTjetMin = min(ptjmin + 5., 1.2 * ptjmin), </i> </li>
<li> <i>MLM:coneRadius = drjmin, </i>
<li> <i>MLM:etaJetMax = etajmax, </i>
</ul>
where the <code>ptjmin</code>, <code>drjmin</code> and
<code>etajmax</code> are the incoming ALPGEN parameters. Note that any
existing values of these parameters are overwritten.
  

<br/><br/><strong>Alpgen:setNjet</strong>  <input type="radio" name="4" value="on" checked="checked"><strong>On</strong>
<input type="radio" name="4" value="off"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>on</strong></code>)<br/>
When switched on, the <code>MLM:nJet</code> parameter (see below) is set
to the incoming <code>njet</code> ALPGEN parameter. Note that any
existing value of this parameter is overwritten.
  


<h2>MLM Merging</h2>

<b>NB: these options are provided by the MLMhooks class,
which must be loaded for this functionality to be present</b>

<p/>
This section describes the MLM merging interface for PYTHIA 8. The most
common reference to the algorithm is [<a href="Bibliography.php" target="page">Man02</a>]. In many respects,
however, the implementation provided in the ALPGEN package should be
considered the official description of the MLM merging procedure.
Although designed primarily to work with events generated with ALPGEN, it
can in principle also be used with events from a different source. This
should not be done without thought, however, and it is up to the user to
understand the details of the algorithm and the implications of using a
different hard process generator.

<p/>
A brief description of the MLM procedure, as implemented, is given here.
First, either the CellJet or SlowJet jet algorithm is chosen, see the
<?php $filepath = $_GET["filepath"];
echo "<a href='EventAnalysis.php?filepath=".$filepath."' target='page'>";?>Event Analysis</a> page for more details.
Both of these algorithms have an <i>R</i> and <i>etaMax</i>
parameter. In addition, CellJet has an <i>eTmin</i> and SlowJet has a
<i>pTmin</i> parameter. These are the primary three parameters of the
merging procedure, and in practice are set dependent on the cuts applied
to the matrix element (ME) generation. We stress that the merging
procedure is not tied to the geometry of a specific physical detector,
but only to the match between the original partons and the resulting
jets, using standard jet algorithms in the phase space region where
partons have been generated.

<p/>
ME samples with different jet multiplicities are run through the event
generator, and the generation interrupted after parton showers have been
applied, but before resonance decays and beam remnants have been
processed. Note in particular that top quarks will not yet
be decayed, which may lead to slight differences with the PYTHIA 6
interface included with the ALPGEN package. In what follows, the
hardness measure of jets/partons is taken to be <i>eT</i> when CellJet
is used and <i>pT</i> when SlowJet is used. The hard system (ignoring
all MPI systems) is then analysed:
<ul>
  <li>
    The particles in the original matrix element process are sorted into
    light partons, heavy partons and other particles. For backwards
    compatibility, a light parton is defined as the set <i>(d, u, s, c,
    b, g)</i> with zero mass. A heavy parton is defined as the set
    <i>(c, b, t)</i> with non-zero mass.
  </li>
  <li>
    All particles not originating from the heavy partons or other
    particles are passed to the jet algorithm and clustered.
  </li>
  <li>
    Clustered jets are matched to the light partons in the original ME
    process. There are two different methods which can be used:
    <ul>
      <li>
        Method 1: The following is done for each parton, in order
        of decreasing hardness. The <i>delta R</i> between the parton
        and all jets is calculated and the smallest value taken. If
        this is less than the jet <i>R</i> parameter, possibly
        multiplied by a constant, the jet and parton are considered to
        match, and the jet is removed from further consideration.
        Note that for CellJet the <i>delta R</i> measure is in
        <i>(eta, phi)</i>, while for SlowJet, it is in
        <i>(y, phi)</i>.
      </li>
      <li>
        Method 2: This method is only possible when using the SlowJet
        algorithm. Before the clustering is performed, extremely soft
        "ghost" particles are added to the event at the
        <i>(y, phi)</i> coordinates of the original matrix element
        partons. If such a particle is clustered into a jet, the parton
        and jet are considered to match. The idea of "ghost" particles
        was originally introduced by FastJet as a way to measure jet
        areas [<a href="Bibliography.php" target="page">Cac06</a>] and should not affect clustering with an
        infrared-safe jet algorithm.
      </li>
    </ul>
  <li>
    If there is a light ME parton remaining which has not been matched
    to a jet, then the event is vetoed. If all ME partons have been
    matched to a jet, but there are still some extra jets remaining,
    then two options are possible:
    <ul>
      <li>
        Exclusive mode: the event is vetoed. This is typically used when
        there are ME samples with higher jet multiplicities, which would
        fill in the extra jets.
      </li>
      <li>
        Inclusive mode: the event is retained if the extra jets are softer
        than the softest matched jet. This is typically used when
        there is no ME sample with higher jet multiplicity, so the parton
        shower should be allowed to give extra jets.
    </ul>
  </li>
  <li>
    All particles originating from the heavy partons are passed to the
    jet algorithm and clustered.
  </li>
  <li>
    The clustered jets are again matched to the original partons, but
    there is no requirement for a match to be present; all matched jets
    are immediately discarded. The matching procedure is much the same
    as for light partons, but with two differences when <i>delta R</i>
    matching is used. First, a different <i>R</i> parameter than that
    used by the jet algorithm may optionally be given. Second, all jets
    that are within the given radius of the parton are matched, not
    just the one with the smallest <i>delta R</i> measure. If there
    are still extra jets remaining then in exclusive mode the event is
    immediately vetoed, while in inclusive mode the event is retained if
    the extra jets are softer than the softest <em>light</em> matched jet.
  </li>
</ul>

<p/>
Some different options are provided, specified further below. These
are set so that, by default, the algorithm closely follows the official
MLM interface provided in the ALPGEN package. All vetoing of events
is done through the usual
<?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?><code>UserHooks</code></a> machinery, and is
therefore already taken into account in the cross section.
In the output from 
<code><?php $filepath = $_GET["filepath"];
echo "<a href='EventStatistics.php?filepath=".$filepath."' target='page'>";?>Pythia::statistics()</a></code>,
the difference between the "Selected" and "Accepted" columns gives the
number of events that have not survived the vetoing procedure. It is
still the responsibility of the user to add together the results from
runs with different jet multiplicities. In the simplest case, when
ALPGEN input is used and the hard process parameters are used to guide
the merging procedure, it is enough to set the <code>MLM:nJetMax</code>
parameter (see below).


<h3>Merging</h3>

<br/><br/><strong>MLM:merge</strong>  <input type="radio" name="5" value="on"><strong>On</strong>
<input type="radio" name="5" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
Master switch to activate MLM merging.
  


<h3>Exclusive mode</h3>

The following settings determine whether clustered jets which do not
match an original hard parton are allowed. They are typically permitted
in the highest jet multiplicity sample, where the parton shower may
produce extra hard jets, without risk of double counting. Any
extra jet produced by the shower must be softer than any matched light
jet, or else the event is vetoed.

<br/><br/><table><tr><td><strong>MLM:exclusive  </td><td>  &nbsp;&nbsp;(<code>default = <strong>2</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Exclusive or inclusive merging.
<br/>
<input type="radio" name="6" value="0"><strong>0 </strong>:  The merging is run in inclusive mode. <br/>
<input type="radio" name="6" value="1"><strong>1 </strong>:  The merging is run in exclusive mode. <br/>
<input type="radio" name="6" value="2" checked="checked"><strong>2 </strong>:  If <ei>MLM:nJet &lt; MLM:nJetMax</ei>, then the merging is run in exclusive mode, otherwise it is run in inclusive mode. If <ei>MLM:nJetMax &lt; 0</ei> or <ei>MLM:nJet &lt; 0</ei>, then the algorithm defaults to exclusive mode. <br/>

<modepick name="MLM:nJet" default="-1" min="-1">
The number of additional light jets in the incoming process.
This is used when <ei>MLM:exclusive = 2</ei>.
Note that for ALPGEN input, this value may be automatically set.

<br/><br/><table><tr><td><strong>MLM:nJetMax  </td><td>  &nbsp;&nbsp;(<code>default = <strong>-1</strong></code>; <code>minimum = -1</code>)</td></tr></table>
This parameter is used to indicate the maximum number of jets that will be
matched. This is used when <ei>MLM:exclusive = 2</ei>.
</modepick>


<h3>Jet algorithm</h3>

The choice of jet algorithm and associated parameters can be adjusted with
the settings below. The PYTHIA 8 internal <code>CellJet</code> and
<code>SlowJet</code> routines are used, see the
<aloc href="EventAnalysis">Event Analysis</aloc> page for more details.

<modepick name="MLM:jetAlgorithm" default="1" min="1" max="2">
The choice of jet algorithm to use when merging against hard partons.
<br/>
<input type="radio" name="7" value="1"><strong>1 </strong>: The CellJet algorithm.<br/>
<input type="radio" name="7" value="2"><strong>2 </strong>: The SlowJet algorithm.<br/>

<br/><br/><table><tr><td><strong>MLM:nEta  </td><td>  &nbsp;&nbsp;(<code>default = <strong>100</strong></code>; <code>minimum = 50</code>)</td></tr></table>
Specific to the CellJet algorithm, the number of bins in pseudorapidity.
</modepick>

<modepick name="MLM:nPhi" default="60" min="30">
Specific to the CellJet algorithm, the number of bins in phi.
</modepick>

<parm name="MLM:eTseed" default="1.5" min="0.0">
Specific to the CellJet algorithm, the minimum <ei>eT</ei> for a cell
to be acceptable as the trial center of a jet.
</parm>

<parm name="MLM:eTthreshold" default="0.1">
Specific to the CellJet algorithm, cells with <ei>eT &lt; eTthreshold</ei>
are completely neglected by the jet algorithm.
</parm>

<modepick name="MLM:slowJetPower" default="-1" min="-1" max="1">
The power to use in the SlowJet algorithm.
<br/>
<input type="radio" name="8" value="-1"><strong>-1 </strong>: The anti-kT algorithm.<br/>
<input type="radio" name="8" value="0"><strong>0 </strong>: The Cambridge/Aachen algorithm.<br/>
<input type="radio" name="8" value="1"><strong>1 </strong>: The kT algorithm.<br/>


<h3>Merging parameters</h3>

The following options are the three main parameters for the merging
procedure. Although here they are in principle free parameters, they should
be heavily influenced by the hard process generation cuts. For ALPGEN
input, these values may be set automatically.

<br/><br/><table><tr><td><strong>MLM:eTjetMin </td><td></td><td> <input type="text" name="9" value="20.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>20.0</strong></code>; <code>minimum = 5.0</code>)</td></tr></table>
For the CellJet algorithm, this gives the minimum transverse energy
inside a cone for a jet to be accepted. For the SlowJet algorithm, this
is instead the minimum transverse momentum required for a cluster to be
accepted as a jet.
  

<br/><br/><table><tr><td><strong>MLM:coneRadius </td><td></td><td> <input type="text" name="10" value="0.7" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>0.7</strong></code>; <code>minimum = 0.1</code>)</td></tr></table>
For the CellJet algorithm, this gives the size of the cone in (eta, phi)
space drawn around the geometric center of the jet. For the SlowJet
algorithm, this gives the <i>R</i> parameter.
  

<br/><br/><table><tr><td><strong>MLM:etaJetMax </td><td></td><td> <input type="text" name="11" value="2.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>2.5</strong></code>; <code>minimum = 0.1</code>)</td></tr></table>
For both jet algorithms, this defines the maximum pseudorapidity that
the detector is assumed to cover. In this context, however, it is tied
to the phase space region in which partons have been generated.
Specifically, particles within <i>etaJetMax + coneRadius</i> are
passed to the jet algorithm, with only jets within <i>etaJetMax</i>
retained in the merging.
  


<h3>Jet matching</h3>

The following parameters control the criteria for matching a clustered jet
to a hard parton.

<br/><br/><table><tr><td><strong>MLM:jetAllow  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 1</code>; <code>maximum = 2</code>)</td></tr></table>
Controls which particles are clustered by the jet algorithm.
<br/>
<input type="radio" name="12" value="1" checked="checked"><strong>1 </strong>:  This option explicitly disallows top quarks, leptons and photons. All other particle types are passed to the jet algorithm. <br/>
<input type="radio" name="12" value="2"><strong>2 </strong>:  No extra particles are disallowed. <br/>

<br/><br/><table><tr><td><strong>MLM:jetMatch  </td><td>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 1</code>; <code>maximum = 2</code>)</td></tr></table>
Criteria for matching a clustered jet to a parton.
<br/>
<input type="radio" name="13" value="1" checked="checked"><strong>1 </strong>:  This option can be used with both the CellJet and SlowJet algorithms. The <ei>delta R</ei> between each parton and jet is taken, and the minimal value compared against <ei>MLM:coneMatchLight * MLM:coneRadius</ei> for light jets or <ei>MLM:coneMatchHeavy * MLM:coneRadiusHeavy</ei> for heavy jets. Note that by default <ei>MLM:coneRadiusHeavy = MLM:coneRadius</ei>, see below. If below this value, the parton and jet are considered to match. With CellJet, the <ei>delta R</ei> measure is in <ei>(eta, phi)</ei>, while with SlowJet it is in <ei>(y, phi)</ei>. <br/>
<input type="radio" name="13" value="2"><strong>2 </strong>:  This option can only be used with the SlowJet algorithm. The hard partons are inserted into the parton level event as "ghost" particles, but at the correct <ei>(y, phi)</ei> position. If this particle is then clustered into a jet, it is considered a match. <br/>

<br/><br/><table><tr><td><strong>MLM:coneMatchLight </td><td></td><td> <input type="text" name="14" value="1.5" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.5</strong></code>; <code>minimum = 0.1</code>)</td></tr></table>
The <code>coneMatchLight</code> parameter used when
<i>MLM:jetMatch = 1</i>.
  

<br/><br/><table><tr><td><strong>MLM:coneRadiusHeavy </td><td></td><td> <input type="text" name="15" value="-1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>-1.0</strong></code>)</td></tr></table>
The <code>coneRadiusHeavy</code> parameter used when
<i>MLM:jetMatch = 1</i>. When assigned a negative value,
the value of <code>MLM:coneRadius</code> is used.
  

<br/><br/><table><tr><td><strong>MLM:coneMatchHeavy </td><td></td><td> <input type="text" name="16" value="1.0" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1.0</strong></code>; <code>minimum = 0.1</code>)</td></tr></table>
The <code>coneMatchHeavy</code> parameter used when
<i>MLM:jetMatch = 1</i>.
  


<h2>Class information</h2>

Some more technical information about the different classes is given
below. For clarity, some limited information on certain private methods
is provided.

<h3>LHAupAlpgen</h3>

This class is derived from the
<?php $filepath = $_GET["filepath"];
echo "<a href='LesHouchesAccord.php?filepath=".$filepath."' target='page'>";?><code>LHAup</code></a> base class, and
uses the standard machinery to pass initialisation and event data to
PYTHIA. These standard functions are not documented here. The complete
parameter file is stored in the PYTHIA Info object, if given, under the
key <code>AlpgenPar</code>.

<a name="method1"></a>
<p/><strong>LHAupAlpgen::LHAupAlpgen(const char *baseFNin, Info *infoPtrIn = NULL) &nbsp;</strong> <br/>
The constructor for the class takes the base filename for the ALPGEN
format files (without file extensions) and optionally a pointer to a
PYTHIA Info class, used for warning/error message printing and for
storing the ALPGEN parameter file. The event and
parameter files are opened immediately, with the <code>AlpgenPar</code>
class, described below, used to parse the parameter file.
  

<a name="method2"></a>
<p/><strong>bool LHAupAlpgen::addResonances() &nbsp;</strong> <br/>
This is a private method used when an event is read in. The information
read from the event file does not always contain a complete listing of
all particles and four-momenta, and so various details must be
reconstructed. Exactly which details are filled in can vary based on the
ALPGEN process in question.
  

<a name="method3"></a>
<p/><strong>bool LHAupAlpgen::rescaleMomenta() &nbsp;</strong> <br/>
This is another private method used when an event is read in.
It shuffles and rescales momenta in an event to ensure energy-momentum
conservation.  First, <i>pT</i> is made to balance by splitting any
imbalance between all outgoing particles with their energies also
scaled. Second, the <i>e/pZ</i> of the two incoming particles are
scaled to balance the outgoing particles. Finally, any intermediate
resonances are recalculated from their decay products.
  

<h3>AlpgenPar</h3>

This class parses an ALPGEN parameter file and makes the information
available through a simple interface. The information is stored
internally in key/value (string/double) format. All lines prior to:
<pre>  ************** run parameters </pre>
are ignored, and in the general case, a line e.g.
<pre>  10   3.00000000000000        ! njets</pre>
would be stored with key "njets" and value "3.0". The following lines
are special cases where the line may be split or the key translated:
<pre>
  3 ! hard process code
  0.000   4.700 174.300  80.419  91.188 120.000 ! mc,mb,mt,mw,mz,mh
  912.905 0.0914176   ! Crosssection +- error (pb)
  100 29787.4  ! unwtd events, lum (pb-1) Njob= 2
</pre>
In the first line, the key "hard process code" is translated to
"hpc". In the second, the mass values are split and each given an entry
in the internal store. In the third, the cross section and cross section
error are stored under the keys "xsecup" and "xerrup" respectively.
Finally, the number of events and luminosity are stored under the keys
"nevent" and "lum" respectively. In the event that a duplicate key is
present, with differing values, the stored value is overwritten and a
warning given.

<a name="method4"></a>
<p/><strong>AlpgenPar::AlpgenPar(Info *infoPtrIn = NULL) &nbsp;</strong> <br/>
The constructor does nothing except for store the PYTHIA Info
pointer, if given. This is used for warning/error message printing.
  

<a name="method5"></a>
<p/><strong>bool AlpgenPar::parse(const string paramStr) &nbsp;</strong> <br/>
This method parses an ALPGEN parameter file. The parameter file is
passed as a single string, mainly intended to be read out from the
PYTHIA Info object using the header information methods.
  

<a name="method6"></a>
<p/><strong>bool AlpgenPar::haveParam(const string &amp;paramIn) &nbsp;</strong> <br/>
Method to check if a parameter with key <code>paramIn</code> is present.
Returns true if present, else false.
  

<a name="method7"></a>
<p/><strong>double AlpgenPar::getParam(const string &amp;paramIn) &nbsp;</strong> <br/>
  
<strong>int AlpgenPar::getParamAsInt(const string &amp;paramIn) &nbsp;</strong> <br/>
Return the parameter with key <code>paramIn</code> as a double or
integer. The presence of a parameter should have already been checked
using the <code>haveParam()</code> function above. If the parameter is
not present, 0 is returned.
  

<a name="method8"></a>
<p/><strong>void AlpgenPar::void printParams() &nbsp;</strong> <br/>
Method to print a list of stored parameters.
  

<h3>AlpgenHooks</h3>

This <?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?><code>UserHooks</code></a> derived class
provides all the <code>Alpgen:*</code> options. It is provided as a
UserHooks class such that the code works regardless of whether ALPGEN
native or LHE file formats are used. It is declared with virtual
inheritance so that it may be combine with other UserHooks classes, see
the "Combining UserHooks" section below.

<a name="method9"></a>
<p/><strong>AlpgenHooks(Pythia &amp;pythia) &nbsp;</strong> <br/>
The constructor takes a PYTHIA object as input, so that the beam
parameter settings can be overridden if the <code>Alpgen:file</code>
option is given. If this is the case, an <code>LHAupAlpgen</code>
instance is automatically created and passed to PYTHIA.
  

<a name="method10"></a>
<p/><strong>bool initAfterBeams() &nbsp;</strong> <br/>
This is the only UserHooks method that is overridden. It is called
directly after PYTHIA has initialised the beams, and therefore the
header information should be present in the PYTHIA Info object. The
<code>AlpgenPar</code> class is used to parse ALPGEN parameters, if
present, which are then used to set further PYTHIA settings.
  
 
<h3>MLMhooks</h3>

This <?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?><code>UserHooks</code></a> derived class
provides all the <code>MLM:*</code> options. It is also declared with
virtual inheritance for combining with other UserHooks classes. It
uses standard UserHooks methods to perform the merging procedure
outlined previously in this manual page, and therefore detailed method
information is not given.

<h3>Combining UserHooks</h3>

It is possible to combine multiple UserHooks classes together, such that
the functionality of both is present. A prerequisite is that the
different UserHooks classes should be declared with virtual inheritance,
e.g.
<pre>
  class MLMhooks : virtual public UserHooks
</pre>
Without this option, when combining two UserHooks derived classes
together, two copies of the base UserHooks class would be created
leading to ambiguity.

<p/>
An example combining the AlpgenHooks and MLMhooks classes together is
given in <code>main32.cc</code>:
<pre>
  class AlpgenAndMLMhooks : public AlpgenHooks, public MLMhooks {
  public:
    AlpgenAndMLMhooks(Pythia &pythia)
      : AlpgenHooks(pythia), MLMhooks() {}
  
    virtual bool initAfterBeams() {
      // Call init of both AlpgenHooks and MLMhooks (order important)
      if (!AlpgenHooks::initAfterBeams()) return false;
      if (!MLMhooks::initAfterBeams())    return false;
      return true;
    }
  };
</pre>
This class inherits from both AlpgenHooks and MLMhooks. Any functions
which are present in both classes should be overridden with a function
that calls the different parent methods in the desired order. In the
above example, the only shared methods are the constructor and
<code>initAfterBeams()</code>.

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

if($_POST["1"] != "void")
{
$data = "Alpgen:file = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "on")
{
$data = "Alpgen:setMasses = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "on")
{
$data = "Alpgen:setMLM = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "on")
{
$data = "Alpgen:setNjet = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "off")
{
$data = "MLM:merge = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "2")
{
$data = "MLM:exclusive = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
if($_POST["7"] != "-1")
{
$data = "MLM:nJetMax = ".$_POST["7"]."\n";
fwrite($handle,$data);
}
if($_POST["8"] != "100")
{
$data = "MLM:nEta = ".$_POST["8"]."\n";
fwrite($handle,$data);
}
if($_POST["9"] != "20.0")
{
$data = "MLM:eTjetMin = ".$_POST["9"]."\n";
fwrite($handle,$data);
}
if($_POST["10"] != "0.7")
{
$data = "MLM:coneRadius = ".$_POST["10"]."\n";
fwrite($handle,$data);
}
if($_POST["11"] != "2.5")
{
$data = "MLM:etaJetMax = ".$_POST["11"]."\n";
fwrite($handle,$data);
}
if($_POST["12"] != "1")
{
$data = "MLM:jetAllow = ".$_POST["12"]."\n";
fwrite($handle,$data);
}
if($_POST["13"] != "1")
{
$data = "MLM:jetMatch = ".$_POST["13"]."\n";
fwrite($handle,$data);
}
if($_POST["14"] != "1.5")
{
$data = "MLM:coneMatchLight = ".$_POST["14"]."\n";
fwrite($handle,$data);
}
if($_POST["15"] != "-1.0")
{
$data = "MLM:coneRadiusHeavy = ".$_POST["15"]."\n";
fwrite($handle,$data);
}
if($_POST["16"] != "1.0")
{
$data = "MLM:coneMatchHeavy = ".$_POST["16"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2012 Torbjorn Sjostrand -->
