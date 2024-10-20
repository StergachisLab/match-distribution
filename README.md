# match-distribution

Sample and force the distribution of PS2 to be the same as that of PS1.

<hr/>
<h4>Software Setup</h4>
Download the software:<br/>
<code>git clone https://github.com/StergachisLab/match-distribution.git</code><br/>
<code>cd match-distribution</code><br/><br/>
Install external dependencies using:<br/>
<code>conda create -n distr-match</code><br/>
<code>mamba env update -n distr-match --file env/matchme.yaml</code><br/><br/>
Activate your environment<br/>
<code>conda activate distr-match</code><br/><br/>
Add needed path to $PATH<br/>
<code>export PATH=${PATH}:`pwd`/src<br/></code>

<hr/>
Usage:<br/>
<code>filter_distr1_distr2.sh model-after.ft.bam input.ft.bam output.filtered.ft.bam</code>

<pre>
model-after.ft.bam:     The distribution you like
input.ft.bam:           The input distribution that will be sampled from and filtered
output.filtered.ft.bam: The output, sampled from Arg2, that has a distribution similar to that of Arg1
</pre>

