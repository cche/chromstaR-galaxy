![chromstaR](https://github.com/ataudt/chromstaR/blob/master/chromstaR_logo.png)
================================

This repository contains a galaxy-wrapper for [chromstaR](https://github.com/ataudt/chromstaR).

Installation
------------

You can add the following lines to your *galaxy/config/tool_conf.xml* file to install it:

   <section name="MyTools" id="mTools">
     <tool file="chromstaR-galaxy/chromstaR_multivariate.xml" />
     <tool file="chromstaR-galaxy/chromstaR_enrichment.xml" />
     <tool file="chromstaR-galaxy/chromstaR_changeCutoff.xml" />
     <tool file="chromstaR-galaxy/chromstaR_differences.xml" />
   </section>


