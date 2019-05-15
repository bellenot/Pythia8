<html>
<head>
<title>The Event Record</title>
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

<form method='post' action='EventRecord.php'>

<h2>The Event Record</h2>

The <code>Event</code> class for event records basically is a vector of 
<code>Particle</code>s, so that it can expand to fit the event size.
The index operator is overloaded, so that <code>event[i]</code> 
corresponds to the <code>i</code>'th particle of an <code>Event</code> 
object <code>event</code>. Thus <code>event[i].id()</code> returns the 
identity of the <code>i</code>'th particle. References to the first, 
<code>i</code>'th and last particle are obtained with 
<code>event.front()</code>, <code>event.at(i)</code> and 
<code>event.back()</code>, respectively.

<p/>
The event size can be found with <code>size()</code>, i.e. valid 
particles are stored in <code>0 &lt;= i &lt; event.size()</code>.

<p/>
Line 0 is used to represent the event as a whole, with its total
four-momentum and invariant mass, but does not form part of the
event history. Lines 1 and 2 contains the two incoming beams, and 
only from here on history tracing works as could be expected. That 
way unassigned mother and daughter indices can be put 0 without 
ambiguity. Depending on the task at hand, a loop may therefore start 
at an index 1 without any loss. Specifically, for translation to other 
event record formats such as HepMC [<a href="Bibliography.php" target="page">Dob01</a>], where the first 
index is 1, the Pythia entry 0 definitely ought to be skipped in order 
to minimize the danger of errors. 

<p/>
New particles can be added to the end of the current event record
with <code>append(Particle)</code>, or 
<code>append(id, status, mother1, mother2, daughter1, daughter2, 
col, acol, p, m, scale)</code>
where <code>p</code> is the four-momentum vector, and everything except 
<code>id</code> defaults to 0. The <code>append</code> method returns 
the index of the new particle position. 
  
<p/>
The existing particle at index <code>iCopy</code> can be copied to the end 
with <code>copy(iCopy, newStatus = 0)</code>. By default, i.e. with 
<code>newStatus = 0</code>, everything is copied precisely as it is, 
which means that history information has to be modified further by hand 
to make sense. With a positive <code>newStatus</code>, the new copy is set 
up to be the daughter of the old, with status code <code>newStatus</code>, 
and the status code of <code>iCopy</code> is negated. With a negative 
<code>newStatus</code>, the new copy is instead set up to be the mother 
of <code>iCopy</code>.

<p/>
A listing of the whole event is obtained with <code>list()</code>. This 
function takes an <code>ostream</code> object as an optional argument.
The basic id, status, mother, daughter, colour, four-momentum
and mass data are always given, but the following switches can be 
set to provide further information, or restrict the set of particles
listed:

<br/><br/><strong>Event:listFinalOnly</strong>  <input type="radio" name="1" value="on"><strong>On</strong>
<input type="radio" name="1" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
The alternative mode, with this option on, omits all intermediate
particles, so that only ones with status > 0 are shown. 
  

<br/><br/><strong>Event:listScaleAndVertex</strong>  <input type="radio" name="2" value="on"><strong>On</strong>
<input type="radio" name="2" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
The alternative mode, with this option on, gives a second line
for each particle, with the production scale (in GeV), the
production vertex (in mm or mm/c) and the invariant lifetime
(also in mm/c).
  

<br/><br/><strong>Event:listMothersAndDaughters</strong>  <input type="radio" name="3" value="on"><strong>On</strong>
<input type="radio" name="3" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
The alternative mode, with this option on, gives a list of all
daughters and mothers of a particle, as defined by the 
<code>motherList(i)</code> and <code>daughterList(i)</code> methods 
described below. It is purely intended for debug purposes, 
e.g. when writing an interface to another event record format. 
  

<br/><br/><strong>Event:extraBlankLine</strong>  <input type="radio" name="4" value="on"><strong>On</strong>
<input type="radio" name="4" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
With this option on, an extra blank line is inserted after each 
particle. There is a tradeoff between improved legibility, in 
particular for <code>Event:listMothersAndDaughters</code> on, and a 
longer listing. 
  

<br/><br/><strong>Event:listJunctions</strong>  <input type="radio" name="5" value="on"><strong>On</strong>
<input type="radio" name="5" value="off" checked="checked"><strong>Off</strong>
 &nbsp;&nbsp;(<code>default = <strong>off</strong></code>)<br/>
With this option on, a list of all junctions in the event is 
displayed after the normal particle list.
  

These choices are stored at the initialization stage of the event.
If you want to change values later, you have to execcute an
<code>Event::initStatic()</code> yourself for them to take effect.  

<p/>
The event record can be emptied for the next event by a 
<code>clear()</code>. The last <code>n</code> entries can be removed by 
<code>popBack(n)</code>, where <code>n = 1</code> by default. 

<p/>
Each particle in the event record has a pointer to the corresponding
particle species in the particle database, used to find some particle
properties. This pointer is automatically set whenever the particle
identity is set by one of the normal methods. Of course its value is 
specific to the memory location of the current run, and so it has no
sense to save it if events are written to file. Should you use some
persistency scheme that bypasses the normal methods when the event is 
read back in, you can use <code>restorePtrs()</code> afterwards to set 
these pointers appropriately.

<p/>
The user would normally be concerned with the <code>Event</code> object 
that is a public member <code>event</code> of the <code>Pythia</code> class. 
Thus, having declared a <code>Pythia</code> object <code>pythia</code>, 
<code>pythia.event[i].id()</code> would be used to return the identity 
of the <code>i</code>'th particle, and <code>pythia.event.size()</code> to 
give the size of the event record. 

<p/>
A <code>Pythia</code> object contains a second event record for the 
hard process alone, called <code>process</code>, used as input for the 
generation of the complete event. Thus one may e.g. call either 
<code>pythia.process.list()</code> or <code>pythia.event.list()</code>. 
To distinguish those two rapidly at visual inspection, the 
<code>"Pythia Event Listing"</code> header is printed out differently, 
in one case adding <code>"(hard process)"</code> and in the other 
<code>"(complete event)"</code>. This is set by a header(string) method.

<p/>
One data member in an <code>Event</code> object is used to keep track of the  
largest <code>col()</code> or <code>acol()</code> tag set so far, so that 
new ones do not clash. The <code>lastcolTag()</code> method returns the 
last tag assigned, i.e. largest value in the current event, and 
<code>nextColTag()</code> ups it by one before returing the value. The 
latter method thus is used when a new colour tag is needed.  

<br/><br/><table><tr><td><strong>Event:startColTag  </td><td></td><td> <input type="text" name="6" value="100" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>100</strong></code>; <code>minimum = 0</code>; <code>maximum = 1000</code>)</td></tr></table>
This sets the initial value used, so that the first one assigned is
<code>startColTag+1</code>, etc. The Les Houches accord [<a href="Bibliography.php" target="page">Boo01</a>] 
suggests this number to be 500, but 100 works equally well.
</modeopen>

<p/>
The <code>scale()</code> methods can be used to set or get the scale 
(in GeV) of the event as a whole. Further methods for event properties 
may be added in the future.

<p/>
There are also a few methods with an individual particle index 
<code>i</code> as input, but requiring some search operations in 
the event record, and therefore not possible to define inside the 
<code>Particle</code> class:

<p/><code>method&nbsp; </code><strong> motherList(i) &nbsp;</strong> <br/>
returns a <code>vector&lt;int></code> containing a list of all the
mothers of a particle. This list is empty for entry 0 , while 
normally it contains one, two or many mothers. The latter case 
applies e.g. to string fragmentation, where the whole fragmenting 
system is counted as mothers. Mothers are listed in ascending order.
  

<p/><code>method&nbsp; </code><strong> daughterList(i) &nbsp;</strong> <br/>
returns a <code>vector&lt;int></code> containing a list of all the
daughters of a particle. This list is empty for a particle that did 
not decay (or, if the evolution is stopped early enough, a parton
that did not branch), while otherwise it can contain a list of 
varying length, from one to many. Many partons may have the same 
<code>daughterList</code>, e.g. in the hard process and fragmentation steps.
For the two incoming beam particles, all shower initiators and beam 
remnants are counted as daughters, with the one in slot 0 being
the one leading up to the hardest interaction. 
  

<p/><code>method&nbsp; </code><strong> iTopCopy(i), iBotCopy(i) &nbsp;</strong> <br/>
are used to trace carbon copies of the particle at index <code>i</code> up 
to its top mother or down to its bottom daughter. If there are no such 
carbon copies, <code>i</code> itself will be returned.
  

<p/><code>method&nbsp; </code><strong> iTopCopyId(i), iBotCopyId(i) &nbsp;</strong> <br/>
also trace top mother and bottom daughter, but do not require carbon 
copies, only that one can find an unbroken chain, of mothers or daughters, 
with the same flavour <code>id</code> code. When it encounters ambiguities,
say a <i>g -> g g</i> branching or a <i>u u -> u u</i> hard scattering,
it will stop the tracing and return the current position. It can be confused
by nontrivial flavour changes, e.g. a hard process <i>u d -> d u</i>  
by <i>W^+-</i> exchange will give the wrong answer. These methods
therefore are of limited use for common particles, in particular for the
gluon, but should work well for "rare" particles. 
  

<p/><code>method&nbsp; </code><strong> sisterList(i) &nbsp;</strong> <br/>
returns a <code>vector&lt;int></code> containing a list of all the 
sisters of a particle, i.e. all the daughters of the first mother,
except the particle itself. 
  

<p/><code>method&nbsp; </code><strong> sisterListTopBot(i) &nbsp;</strong> <br/>
returns a <code>vector&lt;int></code> containing a list of all the 
sisters of a particle, tracking up and back down through carbon copies 
if required. That is, the particle is first traced up with 
<code>iTopCopy()</code> before its mother is found, and then all 
the particles in the <code>daughterList()</code> of this mother are 
traced down with <code>iBotCopy()</code>, omitting the original 
particle itself. Any non-final particles are removed from the list.
Should this make the list empty the search criterion is widened so that
all final daughters are allowed, not only carbon-copy ones. A second
argument <code>false</code> inhibits the second step, and increases 
the risk that an empty list is returned. A typical example of this
is for ISR cascades, e.g. <i>e -> e gamma</i> where the photon 
may not have any obvious sister in the final state if the bottom copy 
of the photon is an electron that annihilates and thus is not part of 
the final state.  
  

<p/><code>method&nbsp; </code><strong> isAncestor(i, iAncestor) &nbsp;</strong> <br/>
traces the particle <code>i</code> upwards through mother, grandmother, 
and so on, until either <code>iAncestor</code> is found or the top of 
the record is reached. Normally one unique mother is required,
as is the case e.g. in decay chains or in parton showers, so that
e.g. the tracing through a hard scattering would not work. For
hadronization, currently first-rank hadrons are identified with the
respective string endpoint quark, which may be useful e.g. for <i>b</i>
physics. However, currently ministrings or junction topologies give 
<code>false</code>. (??)  
  
 
<h3>The Junction Class</h3>

The event record also contains a vector of junctions, which often
is empty or else contains only a very few per event. Methods are
available to add further junctions or query the current junction list.
This is only for the expert user, however, and is not discussed
further here, but only the main points.  

<p/>
A junction stores the properites associated with a baryon number that
is fully resolved, i.e. where three different colour indices are 
involved. There are two main applications,
<ol> 
<li>baryon beams, where at least two valence quarks are kicked out,
and so the motion of the baryon number is notrivial;</li>
<li>baryon-number violating processes, e.g. in SUSY with broken
<i>R</i>-parity.</li>
</ol>
Information on junctions is set, partly in the process generation,
partly in the beam remnants machinery, and used by the fragmentation 
routines, but the normal user does not have to know the details. 

<p/>
For each junction, information is stored on the kind of junction, and 
on the three (anti)colour indices that are involved in the junction.
The possibilities foreseen are:
<ul>
<li><code>kind = 1</code> : incoming colourless particle to three 
outgoing colours (e.g. baryon beam remnant or 
<code>neutralino -> q q q</code>);</li>
<li><code>kind = 2</code> : incoming colourless particle to three 
outgoing anticolours;</li>
<li><code>kind = 3</code> : one incoming anticolor (stored first) 
and two outgoing  colours (e.g. antisquark decaying to quark);</li>
<li><code>kind = 4</code> : one incoming color (stored first) and two 
outgoing anticolours;</li>
<li><code>kind = 5</code> : incoming colour octet to three colours, 
where the incoming colour passes through unchanged and so need not 
be bokkept here, while the incoming anticolor (stored first) and the 
two outgoing colours are (e.g. gluino decay to three quarks);</li>
<li><code>kind = 6</code> : incoming colour octet to three anticolours, 
where the incoming anticolour passes through unchanged and so need not 
be bookkept here, while the incoming color (stored first) and the two 
outgoing colours are.</li>
</ul>
The odd (even) <code>kind</code> codes corresponds to a +1 (-1) change in 
baryon number across the junction.
<br/><b>Warning:</b>Currently only <code>kind = 1, 2</code> are 
implemented.

<p/>
The kind and colour information in the list of junctions can be set 
or read with methods of the <code>Event</code> class, but are not of 
common interest and so not described here.
 
<h3>Subsystems</h3>

The event record also contains a few vectors where parton indices can be 
stored, classified by subsystem. Such information is needed to interleave
multiple interactions, initial-state showers, final-state showers 
and beam remnants. It could also be used in other places. It is intended
to be accessed only by experts, such as implementors of 
<?php $filepath = $_GET["filepath"];
echo "<a href='ImplementNewShowers.php?filepath=".$filepath."' target='page'>";?>new showering models</a>.

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
$data = "Event:listFinalOnly = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "off")
{
$data = "Event:listScaleAndVertex = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
if($_POST["3"] != "off")
{
$data = "Event:listMothersAndDaughters = ".$_POST["3"]."\n";
fwrite($handle,$data);
}
if($_POST["4"] != "off")
{
$data = "Event:extraBlankLine = ".$_POST["4"]."\n";
fwrite($handle,$data);
}
if($_POST["5"] != "off")
{
$data = "Event:listJunctions = ".$_POST["5"]."\n";
fwrite($handle,$data);
}
if($_POST["6"] != "100")
{
$data = "Event:startColTag = ".$_POST["6"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright C 2007 Torbjorn Sjostrand -->
