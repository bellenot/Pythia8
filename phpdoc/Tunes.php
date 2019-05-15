<html>
<head>
<title>Tunes</title>
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

<form method='post' action='Tunes.php'>

<h2>Tunes</h2>

Since some physics aspects cannot be derived from first principles,
this program contains many parameters that represent a true 
uncertainty in our understanding of nature. Particularly afflicted
are the areas of hadronization and multiple interactions, which both
involve nonperturbative QCD physics. 

<p/>
Technically, PYTHIA  parameters can be varied independently of each 
other, but the physical requirement of a sensible description of a set
of data leads to correlations and anticorrelations between the 
parameters. Hence the need to produce tunes, not of one parameter at  
a time, but simultaneously for a group of them. A well-known such 
example is parton densities, where combined tunes to a wide range of data
have been produced, that can then be obtained prepackaged.  

<p/>
Given the many PYTHIA parameters to be tuned, it is convenient to 
divide the task into subtasks. Firstly, if we assume jet universality,
hadronization and final-state parton showers should be tuned to 
<i>e^+e^-</i> annihilation data, notably from LEP1, since this 
offers the cleanest environment. Secondly, with such parameters fixed, 
hadron collider data should be studied to pin down multiple interactions
and other further aspects, such as initial-state radiation. (Thirdly 
would come anything else, such as physics with photon beams, which 
involve further parameters, but that is beyond the current scope.)

<p/>
The first step in this program has now been taken, with a tune to LEP1
data by Hendrik Hoeth, using the Rivet + Professor framework. Starting
with version 8.125 it defines the default values for hadronization 
parameters and timelike showers. 

<p/>
The situation is worse for multiple interactions, where PYTHIA 8 is more 
different from PYTHIA 6. Nevertheless, a first simple tune is now 
available, appropriately called "Tune 1", and is default starting with 
version 8.127.

<p/>
In the future we hope to see further PYTHIA 8 tunes appear. Like with 
parton distributions, there is likely to be several tunes, because 
different sets of data will pull in different directions, by imperfections   
in the model or in the data, and by differences in the chosen
tuning strategies. We therefore propose to collect some of these tunes
here, in a prepackaged form. Of course, in all cases it is a matter
of setting values for parameters already defined elsewhere, so the
tunes offer no new functionality, only a more convenient setup. 

<p/>
If you set either the <code>Tune:ee</code> and <code>Tune:pp</code> 
modes below non-zero then all parameters used in the respective tune 
will be set accordingly when <code>pythia.init(...)</code> is called. 
You can check this by calling <code>pythia.settings.listChanged()</code> 
before and after initialization; before only the tune modes are 
nondefault, afterwards all the non-default-valued parameters in the 
tune appear. Therefore, for better or worse, you cannot combine a tune 
option with your own choices for some of the parameters used in the tune, 
since the values you set before <code>pythia.init(...)</code> would be 
overwritten at that point. 

<br/><br/><table><tr><td><strong>Tune:ee  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 3</code>)</td></tr></table>
Choice of tune to <ei>e^+e^-</ei> data, mainly for the hadronization
and timelike-showering aspects of PYTHIA. 
<br/>
<input type="radio" name="1" value="0" checked="checked"><strong>0 </strong>: no values are overwritten at initialization,  so you can set the individual parameters as you wish. <br/>
<input type="radio" name="1" value="1"><strong>1 </strong>: the original PYTHIA 8 parameter set, based on some very old flavour studies (with JETSET around 1990) and a simple tune  <ei>of alpha_strong</ei> to three-jet shapes to the new  <ei>pT</ei>-ordered shower. These were the default values before version 8.125.  <br/>
<input type="radio" name="1" value="2"><strong>2 </strong>: a tune by Marc Montull to the LEP 1 particle composition, as published in the RPP (August 2007). No related (re)tune  to event shapes has been performed, however.   <br/>
<input type="radio" name="1" value="3"><strong>3 </strong>: a tune to a wide selection of LEP1 data by Hendrik  Hoeth within the Rivet + Professor framework, both to hadronization and timelike-shower parameters (June 2009). These are the default values  starting from version 8.125, so currently there is no need for this option. <br/>

<br/><br/><table><tr><td><strong>Tune:pp  </td><td>  &nbsp;&nbsp;(<code>default = <strong>0</strong></code>; <code>minimum = 0</code>; <code>maximum = 2</code>)</td></tr></table>
Choice of tune to <ei>pp / ppbar</ei> data, mainly for the 
initial-state-radiation, multiple-interactions and  beam-remnants
aspects of PYTHIA. 
<br/>
<input type="radio" name="2" value="0" checked="checked"><strong>0 </strong>: no values are overwritten at initialization,  so you can set the individual parameters as you wish. <br/>
<input type="radio" name="2" value="1"><strong>1 </strong>: old default used up to version 8.126, based on  some early and primitive comparisons with data. <br/>
<input type="radio" name="2" value="2"><strong>2 </strong>: "Tune 1", new default from 8.127 onwards,  based on some data comparisions by Peter Skands. <br/>



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

if($_POST["1"] != "0")
{
$data = "Tune:ee = ".$_POST["1"]."\n";
fwrite($handle,$data);
}
if($_POST["2"] != "0")
{
$data = "Tune:pp = ".$_POST["2"]."\n";
fwrite($handle,$data);
}
fclose($handle);
}

?>
</body>
</html>

<!-- Copyright (C) 2010 Torbjorn Sjostrand -->
