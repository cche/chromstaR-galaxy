![chromstaR](https://github.com/ataudt/chromstaR/blob/master/chromstaR_logo.png)
================================

This repository contains a galaxy-wrapper for [chromstaR](https://github.com/ataudt/chromstaR).

Installation
------------

Clone this repository into *galaxy/tools/* with

```
git clone https://github.com/cche/chromstaR-galaxy.git
```

Then add the following lines to your *galaxy/config/tool_conf.xml* file to install it:

```xml
<section name="MyTools" id="mTools">
  <tool file="chromstaR-galaxy/chromstaR_multivariate.xml" />
  <tool file="chromstaR-galaxy/chromstaR_enrichment.xml" />
  <tool file="chromstaR-galaxy/chromstaR_changeCutoff.xml" />
  <tool file="chromstaR-galaxy/chromstaR_differences.xml" />
</section>
```

