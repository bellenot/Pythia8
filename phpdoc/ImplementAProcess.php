<html>
<head>
<title>Implement a Process</title>
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

<form method='post' action='ImplementAProcess.php'>

<h2>Implement a Process</h2>

Normally users should not try to implement new processes internally
inside <code>Pythia</code>, but rather use the 
<?php $filepath = $_GET["filepath"];
echo "<a href='LesHouchesAccord.php?filepath=".$filepath."' target='page'>";?>Les Houches Accord</a> 
to input your own events for processing by <code>Pythia</code>. 
However, in case you do want to implement a new process at the same 
level as the internal <code>Pythia</code> ones, here is a brief summary. 
For the details, of course, it always makes sense to take a similar
already-implemented process as a starting point.

<p/>
There are two steps involved in implementing a process:
<br/>1) implementing the matrix elements, including information on
incoming and outgoing flavours and colours, and 
<br/>2) making the process available to the rest of the program,
so that it can be switched on by a user.
<br/>We consider these two aspects in turn.

<h3>Matrix Elements</h3> 

The matrix-element information have to be stored in a matched pair
of <code>.cc</code> and <code>.h</code> files, in the 
<code>src</code> and <code>include</code> directories, respectively.
These could be one set of existing files (<code>SigmaQCD</code>,
<code>SigmaEW</code>, <code>SigmaOnia</code> and 
<code>SigmaSUSY</code>), if it belongs there, or entirely new files, 
say <code>SigmaCompositeness</code>, where the name attaches to the 
scenario studied.

<h4>Header</h4>

<p/>
The <code>.h</code> file is used to define the class for the new 
process. This class has to be derived either from  
<code>Sigma1Process</code>, for <i>2 -> 1</i> processes, or from 
<code>Sigma2Process</code>, for <i>2 -> 2</i> ones. (The 
<code>Sigma0Process</code> class is used for elastic, diffractive
and minimum-bias events, and should not be used.) Both
<code>Sigma1Process</code> and <code>Sigma2Process</code>
are in their turn derived from the <code>SigmaProcess</code> 
base class.

<p/>
The header information of a new class must contain a constructor. 
The constructor can take arguments to allow a set
of related processes to share common code. For instance,
g g -> Q Qbar, Q = c or b, is only coded once, and then the 
constructor takes the quark code (4 or 5)  as argument, 
to allow the proper amount of differentiation. A destructor is only
needed if you use <code>new</code> yourself. (The <code>InFlux</code> 
object, see below, is automatically taken care of.)


<p/>
A few simple methods are used to encapsulate information on the 
particular process:
<br/>* <code>name()</code> returns the name of the process, as it will 
be shown in listings.
<br/>* <code>code()</code> returns an integer identifier of the process.
This has no internal function, but is only intended as a service for
the user to rapidly identify which process occured in a given event.
<br/>* <code>id3Mass(), id4Mass()</code> are the one or two final-state 
flavours, where masses are to be selected before the matrix elements
are evaluated. Only the absolute value need to be given. For massless
particles, like gluons and photons, one need not give anything, i.e. 
one defaults to 0. The same goes for normal light quarks, where masses 
presumably are not implemented in the matrix elements.  Later on, these 
quarks can still (automatically) obtain constituent masses, once a 
u, d or s flavour has been selected. 
<br/>* <code>resonanceA(), resonanceB()</code> are the codes of up to two 
s-channel resonances contributing to the matrix elements. These are 
used by the program to improve the phase-space selection efficiency,
by partly sampling according to the relevant Breit-Wigner. Massless
resonances need not be specified.

<p/>
In addition, the <code>initProc()</code>, <code>initFlux()</code>, 
<code>sigmaHat()</code> and <code>setIdColAcol()</code> methods have to be 
implemented. This often requires more code, and if so the implementation 
belongs in the <code>.cc</code> file.

<h4>Initialization</h4>

<code>initProc()</code> and <code>initFlux()</code> are called
during initialization and then not again. The processes used in the 
multiple-interactions framework do not call <code>initFlux()</code> 
since the parton densities are handled centrally in this case.  

<p/>
<code>initProc()</code> is the place where all types of process-specific 
initialization calculations could be put, e.g. readout and storage of
masses and couplings. Many processes have no such object.

<p/>
<code>initFlux()</code> is not required for the "inclusive" processes,
i.e. the elastic, diffractive and minimum-bias ones. In all other cases  
this object has to specify the combinations of incoming partons that are 
allowed for the process under consideration. This is obtained by creating 
a pointer to an <code>InFlux</code> object. Currently allowed options are 
<br/>* <code>InFluxgg</code>: g g , 
<br/>* <code>InFluxqg</code>: q g or qbar g,
<br/>* <code>InFluxqqbarqqDiff</code>: q q', q qbar' or qbar qbar', 
q and q' different flavours,
<br/>* <code>InFluxqqDiff</code>: q q' or qbar qbar', 
q and q' different flavours, 
<br/>* <code>InFluxqqSame</code>: q q or qbar qbar, same flavour,
<br/>* <code>InFluxqqbarDiff</code>: q qbar', q and q' different 
flavours,
<br/>* <code>InFluxqqbarSame</code>: q qbar, same flavour,
<br/>* <code>InFluxff</code>: f f', f fbar' or fbar fbar', where
f and f' may be same or different,
<br/>* <code>InFluxffbarSame</code>: f fbar, same flavour,
<br/>* <code>InFluxffbarChg</code>: f fbar', where the combination 
has charge +-1 (e.g. u dbar, e+ nu_e).
<br/>For a hadron beam, the generic fermion alternatives default 
back to quarks only.

<p/>
The <code>InFlux</code> object is created with <code>new</code> but
no special desctructor is needed, since this is handled by the 
base class destructor.

<p/>
After the allowed channels have been set up, it is possible to
specify that the channels come with a different weight, apart from
the automatic convolution with the respective parton densities
<br/>* <code>weightCharge2()</code>: weights by the squared quark 
charge (also works e.g. for q g states, but not for a quark pair 
of different charge).
<br/>* <code>weightCKM2()</code>: weights by the squared of the 
CKM mixing-matrix element, vanishing when no such element exists. 
<br/>* <code>weightCKM2sum(int mode, int idQ)</code>: provides
the CKM-related weight for a number of different situations.
<br/>1) Sum of CKM weights for transformation of an incoming flavour
to an outgoing one, excluding top (used e.g. for f g -> f' W).
<br/>2) Ditto, but product of both incoming sides (not used currently).
<br/>3) Coupled CKM weights of the two sides, as consistent with
t-channel W exchange (used e.g. for q_1 q_2 -> q_3 q_4).
<br/>4) Coupled CKM weights of the two sides, but on one side to
specified flavour idQ (used e.g. for q_1 q_2 -> t q_3).
<br/>* <code>weightInvCol()</code>: factor 1/3 for an incoming
qqbar pair and 1/8 for an incoming g g one.
<br/>* <code>weightNeutrinoSpin()</code>: multiply by factor 2 
for an incoming neutrino, since they always are lefthanded and so
are not averaged over incoming spin.
<br/>* <code>weightFixed(double nowWeight)</code>:
weights all channels by <code>nowWeight</code>.
<br/>* <code>weightFixed(int id1, int id2, double nowWeight, 
bool flipSide = true, bool conjugate = true, bool allGen = true)</code>:
weights by <code>nowWeight</code> for the specific incoming state
consisting of  <code>id1</code> and <code>id2</code>. If
<code>flipside</code> then also apply this weight when the two
incoming sides are interchanged (but only once for 
<code>id1 = id2</code>, of course). If <code>conjugate</code>
then also apply it when quarks and leptons are replaced by their
antiparticles. If <code>allGen</code> then apply weights provided
for first-generation fermions also to equivalents in the second
and third generation (i.e. u dbar is also applied to u sbar, u bbar,
c dbar, c sbar and c bbar).      

<p/>
You can see the list of allowed channels and their respective fixed 
weights by switching on the <code>InFlux:showChannels</code> flag.
Note that channels which are assigned a vanishing weight during the
initialization step are pruned from the channel list. 

<h4>Event weight</h4>

<code>sigmaHat()</code> encodes the relevant matrix element. For a 
<i>2 -> 1</i> process, what is to be provided is 
<i>sigmaHat(sHat)</i>, where <code>sH</code> and <code>sH2</code>
are available to be used. For a <i>2 -> 2</i> process, instead
<i>d(sigmaHat)/d(tHat)</i> should be calculated, based on
provided  <code>sH, sH2, tH, tH2, uH, uH2, m3, s3, m4, s4</code> and 
<code>pT2</code> values (<code>s3 = m3*m3</code>). In either case,
<i>alpha_s</i> and <i>alpha_em</i> have also been calculated, 
and are stored in <code>alpS</code> and <code>alpEM</code>. Also 
other standard variables may be used, like 
<code>CoupEW::sin2thetaW()</code>, and related flavour-dependent
vector and axial couplings in <code>CoupEW</code> and CKM combinations
in <code>VCKM</code>. Note that, normally the cross sections would come 
in dimensions of inverse GeV to second or fourth power, respectively. 
To translate it to the compulsory mb or mb/GeV^2 dimensions, you need 
to multiply by the <code>CONVERT2MB</code> constant.

<p/>
A further machinery exists when event-by-event channel-dependent 
weights are required. A typical example would be that,
for q qbar -> gamma*/Z0, the relative weight of incoming u-type and 
d-type quarks depends on the relative admixture of gamma* and of Z0, 
which depends on the subprocess energy, and so varies from one event 
to the next. Therefore, in addition to the fixed weight above, to be 
set at initialization, each channel contains a variable weight factor, 
and the total weight of a channel is a product of the two. The variable 
weight is set to 1 at the creation of a new channel, so if not touched 
it makes no difference. It is not changed between events, however, 
i.e. not automatically reset to 1. It can be manipulated by two methods: 
<br/>* <code>weightInState(double nowWeight = 1.)</code>: set all variable 
weights to the value specified (usually not necessary).
<br/>* <code>weightInState(int id1, int id2, double nowWeight, 
bool flipSide = true, bool conjugate = true, bool allGen = true)</code>:
set the variable weight factor for the incoming state defined by 
<code>id1</code> and <code>id2</code> to be <code>varWeight</code>
Se above for an explanation of the <code>flipside</code>,
<code>conjugate</code> and <code>allGen</code> options.  

<p/>
Normally, only one cross section is returned from <code>sigmaHat()</code>,
and the standard machinery takes care of weighting it with parton
densities and the other factors specified above. However, if you use 
variable weights, you should use the <code>weightInState(...)</code> 
method above to specify the weight for each incoming flavour combination.
Note that, when weights are used, there is a complete freedom how you 
shuffle contributions between the individual channels and the overall 
weight returned by <code>sigmaHat()</code>, since only the product is 
relevant. For sanity, as much as possible should be put in the 
<code>sigmaHat()</code>, with <code>weightFixed(...)</code> and
<code>weightInState(...)</code> only providing dimensionless 
extrafactors of order unity.

<p/>
Processes that should also be used by the multiple-interactions 
machinery have to be coded somewhat differently in cases where 
the weight depends on the incoming flavour combination. Since no 
<code>InFlux</code> object is created in this case it does not help
to impose channel-by-channel preweights. On the other hand, when a
process is called from <code>MultipleInteractions</code> the two 
incoming flavours are already known, so there is no need to loop over
flavour combinations with different weights. In a run a given process
can be used both as a hard process and as part of the multiple 
interactions, but these are created as two separate objects. These
are distinguished by the variable <code>id12IsSet = false</code> in
the former case and  <code>= true</code> in the latter, so some small
piece of <code>sigmaHat()</code> may need to distinguish these two
possibilities. An example would be <i> q qbar -> g gamma</i>,
where the weighting by squared charge could be taken care of once and
for all in <code>initFlux</code> for a hard process, whereas it would
have to be included explicitly in <code>sigmaHat()</code> for a 
multiple-interactions process, making use of the fact that <code>id1</code> 
(and <code>id2 = -id1</code>) are already set in that case. However, 
users are strongly urged not to try to include processes of their own
among multiple interactions, so this paragraph is mainly to explain
the meaning of some extra lines of code for some existing processes.

<h4>Finalize event</h4>

<code>setIdColAcol()</code> is called when a kinematical configuration
has been accepted and flavours and colours need be provided. The first
part is handled by the <code>setId</code> method. When this routine is 
called the two incoming flavours have already been set in <code>id1</code> 
and <code>id2</code>, whereas the one or two outgoing ones either are fixed 
for a given process or can be determined from the instate (e.g. whether a 
<i>W+</i> or <i>W-</i> was produced). The colours are handled by the 
<code>setColAcol</code> method. Les Houches style colour tags are used, 
but starting with number 1. The input is grouped particle by particle, 
with the colour index before the anticolour one. You may need to select 
colour flow dynamically, depending on the kinematics, when several distinct 
possibilities exist. Trivial operations, like swapping colours and
anticolours, can be done with existing methods. There is also a standard 
method to pick a final flavour from an initial one with CKM mixing.

<p/>
When the <code>id3Mass()</code> and <code>id4Mass()</code>
methods have been used, the order of the outgoing particles may be 
inconsistent with the way the <i>tHat</i> and <i>uHat</i>
variables have been defined. A typical example would be a process like
<i>q g -> q' W</i> with <i>tHat</i> defined between incoming and 
outgoing quark, but where <code>id3Mass() = 24</code> and so the
process is to be stored as <i>q g -> W q'</i>. One should then put
the variable <code>swapTU = true</code> for each event where the 
<i>tHat</i> and <i>uHat</i> variables should be swapped before 
the event kinematics is reconstructed. This variable is automatically
restored to <code>false</code> for each new event.

<h3>Access</h3> 

A flag has to be defined, that allows the process to be switched on;
by default it should always be off. The name of the flag should be 
chosen of the type <code>model:process</code>. Here the 
<code>model</code> would be related to the general scenario considered,
e.g. <code>Compositeness</code>, while <code>process</code> would
specify instate and outstate, separated by a 2 (= to), e.g.
<code>ug2u*g</code>. 
When several processes are implemented and "belong together" it is 
also useful to define a <code>model:all</code> switch that affects
all the separate processes. 

<p/>
The flags should normally be stored in the <code>Processes.xml</code>
file. This is to make them easily found by users. You could create
and use your own <code>.xml</code> file, so long as you then add that 
name to the list of files in the <code>Index.xml</code> file. (If not,
the flags would never be created and the program would not work.)  

<p/>
In the <code>ProcessContainer.c</code> file, the final 
<code>SetupContainers::init</code> routine needs to be expanded to
create instances of the processes switched on. This code is fairly
repetitive, and should be easy to copy and modify from the code 
already there. The basic structure is 
<br/>(i) check whether a process is requested by the user and, if so, 
<br/>(ii) create an instance of the matrix-element class, 
<br/>(iii)create a container for the matrix element and its associated 
phase-space handling, and 
<br>(iv) add the container to the existing process list.  

<p/>
Two minor variations are possible. One is that a set of related 
processes are lumped inside the the same initial check, i.e. are 
switched on all together. The second is that the matrix-element 
constructor may take arguments, as specified by you (see above). 
If so, the same basic matrix element may be recycled for a set of 
related processes, e.g. one for a composite u and one for a composite d. 
Obviously these variations may be combined.

</body>
</html>

<!-- Copyright C 2007 Torbjorn Sjostrand -->
