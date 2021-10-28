############################################################
# For personal site, 03.20
# This generates the file "index.html"
############################################################

import sys, os
sys.path.append('..')
import lib.read_chunks as RC

######################
# HTML template
######################

html_template = """
<!doctype html>
    {head}

<body>
	<div class="row" id="top_grid">
		<div class="col-2-24" id="margin"></div>
		<div class="col-20-24" id="main_header">Read mapping simulations</div>
		<div class="col-2-24" id="margin"></div>
	</div>

	{nav}


	<div class="row">
		<div class="col-2-24" id="margin"></div>
		<div class="col-20-24" id="main_col">
			<h2>
				Useful software
			</h2>

            <ul>
				<h3>
					<a href="https://github.com/ncsa/NEAT" target="_blank">NEAT read simulator</a>
				</h3>
			<ul>
		</div>
        <div class="col-2-24" id="margin"></div>
	</div>

	<div class="sep_div"></div>

	<div class="sep_div"></div>
	<div class="sep_div"></div>
	<div class="sep_div"></div>
	<div class="sep_div"></div>
	<div class="sep_div"></div>
	<div class="sep_div"></div>
	<div class="sep_div"></div>
	<div class="sep_div"></div>

    {footer}
</body>
"""

######################
# Main block
######################
pagefile = "links.html";
print("Generating " + pagefile + "...");
title = "MapSims - links"

head = RC.readHead(title, pagefile);
nav = RC.readNav(pagefile);
footer = RC.readFooter();

outfilename = "../../" + pagefile;

with open(outfilename, "w") as outfile:
    outfile.write(html_template.format(head=head, nav=nav, footer=footer));