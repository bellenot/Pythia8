<html>
<head>
<title>Event Statistics</title>
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

<form method='post' action='EventStatistics.php'>

<h2>Event Statistics</h2>

At the end of the run you will want to write out the final statistics
on number of events generated, the corresponding cross sections and 
the number of errors encountered. This is done with
<pre>
    pythia.statistics();
 </pre>
assuming <code>pythia</code> is an instance of the <code>Pythia</code>
class. The <code>statistics</code> method in its turn calls on the 
methods below.

<h3>Cross-section statistics</h3>

The <code>ProcessLevel::statistics()</code> member will loop over the
list of existing processes, and for each write out name, code,
the number of tried, selected and accepted events, the cross section and 
the estimated error on the latter. The three different event numbers are 
related to the Monte Carlo method used, whereby an initial upper estimate
of the cross section is used to select a large number of trial phase-space 
points, whereof then not all survive. Rejections are normally done by the
internal machinery, but can also be obtained by
<?php $filepath = $_GET["filepath"];
echo "<a href='UserHooks.php?filepath=".$filepath."' target='page'>";?>user hooks</a>. 
Therefore:<br/>
(i) tried events reflect the original number of phase-space points 
probed, as part of the upper estimate;<br/>
(ii) selected events correspond to those that survive the internal 
Monte-Carlo selection procedure;<br/> 
(iii) accepted events are those that also survive the additional 
user cuts.<br/> 
In most runs there would be no user hooks implemented, and then the 
numbers of selected and of accepted events will agree.

<h3>Error messages</h3>

When Pythia is run, errors may occur, and give rise to warning messages.
These may be of varying severity, as follows:
<br/><b>Abort</b> means things went seriously wrong, and the 
initialization or event generation failed. In the former case it is 
not possible to generate events at all, in the latter the current
event is flawed and should be skipped. In either case the respective
method, <code>pythia.init(...)</code> or <code>pythia.next()</code>,
then also returns the value <code>false</code>. There are occasions
where an abort may be deliberate, such as when a file of Les Houches
Events is read and the end of the file is reached.
<br/><b>Error</b> normally is less severe. Typically the program will
back up one step and try again. There are cases where this is not possible,
in particular during the initialization and the generation of a hard
process, and then the error may be followed by an abort as a direct 
consequence (with two separate messages).   
<br/><b>Warning</b> is even less severe. In some cases the program will 
try again, with  good chances of success, in others no measure at all
need to be taken. 

<p/>
The <code>ErrorMsg</code> class is rather small. It is handed any 
abort, error or warning messages during the event generation phase, and 
will store each distinct message, with a counter for how many times it is 
issued. Thus it is possible to limit the number of identical messages 
issued. The summary table printed by <code>pythia.statistics()</code> 
provides a table with all the different messages issued, in 
alphabetical order, with the total number of times each was generated.

<p/>
Should you use several <code>pythia.init(...)</code> calls in your code,
the error counters will be reset for each new one. 

<p/>
There is only one mode affecting the <code>ErrorMsg</code> operation:
<br/><br/><table><tr><td><strong>ErrorMsg:timesToPrint  </td><td></td><td> <input type="text" name="1" value="1" size="20"/>  &nbsp;&nbsp;(<code>default = <strong>1</strong></code>; <code>minimum = 0</code>)</td></tr></table>
The number of times each distinct message is printed. That is, by
default, each new kind of error/warning is only printed once.
</modeopen> 

<h3>Multiple-interactions statistics</h3>

If you call <code>pythia.statistics(true)</code>, i.e. with the optional
argument <code>true</code>, also statistics on multiple interactions
is printed, comprising a list of all allowed subprocesses with how
many times each of them has been generated. For the minimum-bias
process this also includes the hardest interaction, while else the 
hardest process is excluded from the statistics. (This is because 
the hardest process is of the same character and generated by the same
machinery in the former case but not in the latter. Also, for the 
former case only, the standard statistics listing only lists 
minimum bias as one single process, i.e. does not further specify 
the character of the hardest subprocess, so there is not any overlap 
between the two.)

<p/>
Should you use several <code>pythia.init(...)</code> calls in your code,
the multiple-interactions statistics will be reset for each new one. 

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

if($_POST["1"] != "1")
{
$data = "ErrorMsg:timesToPrint = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2007 Torbjorn Sjostrand -->
