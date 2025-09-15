# qp-template
Project directory and LaTeX template for UBC STAT 548 Qualifying Paper reports

## Instructions for use:

1. Click "Use Template" at the top. This basically forks it into your own GitHub
   account. Don't have a GitHub, set one up now <https://github.com>. Note: this
   is hosted on servers outside Canada.
1. If necessary, [install git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).
1. Clone your new repo to your machine somewhere.
1. At some point before you submit the report,
   [add me as a collaborator](https://help.github.com/en/articles/inviting-collaborators-to-a-personal-repository).<sup id="a1">[1](#f1)</sup>
1. Use the directory structure to keep your project organized and use the LaTeX
   template in the `report` directory for your report. See below for more
   details on using the LaTeX template. See 
   [my website](https://dajmcdon.github.io/teaching/stat548.html) 
   for details on my expectations for the report.
1. Commit and push your changes regularly. It serves as a back-up and lets me
   see how you progressed.
1. When you are ready to submit your report, commit with the message "submission
   for DJM review".
1. Send me an email confirming that you have submitted your report for review
   and marking.


## Details on directory usage

The short version: Each subdirectory has its own README file that tells you what goes in the directory.

The long(er) version: 

* Use the `report/` directory to for the written portion of the project. 
* Use `report/qp.bib` to keep your biliography up-to-date, and use it as the
  bib-file for your report. 
* Use `src/` for any code you develop and `data/` for any data you use
  (including output from your code). (There may be good reasons for not keeping
  your code in the same repository as your report, in which case let me know and
  invite me as a collaborator to that repo, as well.)

## Details on LaTeX template usage

Please use the LaTeX template in the `report/` directory for your report. It is
set up to use the `rss` bibliography style, which is what I prefer. You can
compile the report using `pdflatex` and `bibtex` from the command line, or your
preferred LaTeX editor.

I personally avoid Overleaf because it costs money to use useful features, and
doesn't operate as smoothly as Emacs or VSCode.

Please read through the comments in the LaTeX file for more details on how to
use it.

_Acknowledgements_: [Dave Blei](http://www.cs.columbia.edu/~blei/) for (most of)
the basic directory structure and reading/project logs. As should be clear from
the fork, I stole this from [Ben](https://www.stat.ubc.ca/~benbr/)

<b id="f1">[1]</b> *When* you give me access is entirely up to you. I will not
look at your repository before your submit your report unless you ask me to (and
even then I may not). [↩](#a1)
