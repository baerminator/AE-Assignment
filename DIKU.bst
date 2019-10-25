% Release $Date: 2018/11/06 17:19:53 $
%
% Comments can be sent to jyrki@di.ku.dk
%
% Development history
%
% 2018-11-02 Unified url support -Jyrki
%   Added doi support
%   Added isbn support
%   Added eprint support
%   Added fields: doi, isbn, isbn-13, url, lastaccessed, articleno, numpages, 
%   archiveprefix, eprint, eprinttype, eprintclass, primaryclass
%
% 2008-05-15 Now techrept accepts an editor -Jyrki
% 2007-10-21 Case change in the thesis type omitted -Jyrki
% 2005-02-04 Now author names are output in roman, not in small caps -Jyrki
% 2004-01-27 Organization added to book -Jyrki
% 2002-10-10 Now url is printed by techreport. -S�ren & Jyrki
% 2002-08-29 Now unpublished outputs the title in roman. -Jyrki
% 2002-08-15 I corrected an error in the formatting of "von" names. -Jyrki
% 2001-10-12 Jeppe wanted that editor is removed inproceedings. He asked whether
%            url should be provided by other entries than webpages. -Jyrki
% 2001-04-17 article, techreport and thesis titles in lower-case letters -Jyrki
% 2000-09-26 Use url package to render url's -Jyrki
% 2000-07-03 Testing and correcting some small errors -Jyrki
% 2000-06-30 Minor revisions in series and inproceedings -Jyrki
% 2000-06-15 Changes in bvolume -Jyrki
% 2000-06-07 Changes in article -Jyrki
%
% Adaption by Tomi Pasanen
%	of
% ACM Transactions bibliography style (24-Jan-88 version)
% filename "acm.bst"
%	and of
% American Mathematical Society bibliography style 13-Oct-1994 version 1.2beta
% filename "amsplain.bst"
%
% Major rewrites by Jyrki
% - 2000-06-30
% - 2018-11-02
%
% Status: bibliography entries 
% booklet,
% inbook,
% misc,
% have not been tested.
% Cross-references have never been used.
%
% BibtTeX is a primitive stack engine created by Oren Patashnik. A style 
% file like this translates the bibliography entries mentioned in an aux-file 
% and described in a bib-file to a sequence of tokens (name, title, journal, 
% etc.) in the LaTeX-format to a bbl-file. This program considers the entries
% one at a time, generates the tokens for each entry, and the mode tells
% which separator is outputted _before_ each token. The default separator is
% blank (" "); other separators are comma (", "), colon (": "), and period
% (". "). Knowing this it should be relatively easy to modify the output
% corresponding to any bibliography entry if necessary.
%
% This style file also knows non-standard entries webpage and website
% which have the field url giving the location of the resource.
%
% Accepted by all entries: doi, eprint, eprinttype, eprintclass, primaryclass, 
% archiveprefix, isbn, isbn-13, url, lastaccessed
%
% article
% Required: author, title, journal, year, volume
% Optional: number, pages, month, doi, note, key
% Non-standard: articleno, numpages; month not rendered
%
% book
% Required: author/editor, title, publisher, year
% Optional: volume/number, series, address, edition, month, note, key, url
% Non-standard: organization, isbn, isbn-13; address & month not rendered
%
% booklet
% Required: title
% Optional: author, howpublished, address, month, year, note, key
% Non-standard: month not rendered
%
% inbook
% Required: author/editor, title, chapter/pages, publisher, year
% Optional: volume/number, series, type, address, edition, month, note, key
% Non-standard: month not rendered
%
% incollection
% Required: author, title, booktitle, publisher, year
% Optional: editor, volume/number, series, type, chapter, pages, address, 
% edition, month, note, key
% Non-standard: month not rendered
%
% inproceedings
% Required: author, title, booktitle, year
% Optional: editor, volume/number, series, pages, address, month, organization,
% publisher, note, key
% Non-standard: editor & month not rendered; address refers to the address of 
% the publisher, not to the location where the confrerence was held
%
% mastersthesis
% Required: author, title, school, year
% Optional: type, address, month, note, key
% Non-standard: month not rendered
%
% phdthesis
% Required: author, title, school, year
% Optional: type, address, month, note, key
% Non-standard: month not rendered
%
% misc
% Required: none
% Optional: author, title, howpublished, month, year, note, key
% Non-standard: month not rendered
%
% techreport
% Required: author, title, institution, year
% Optional: type, number, address, month, note, key
% Non-standard: editor, eprint, eprinttype, eprintclass, primaryclass;
% month not rendered
%
% unpublished
% Required: author, title, note
% Optional: month, year, key
% Non-standard: month not rendered
%
% webpage
% Required: (author or editor) or (organization and key), title, url, year
% Optional: type, note
%
% website
% Required: editor or organization (or both), key, title, url, year.
% Optional: note

ENTRY
  { address
    archieveprefix
    articleno
    author
    booktitle
    chapter
    edition
    editor
    howpublished
    institution
    journal
    key
    month
    note
    number
    numpages
    organization
    pages
    publisher
    school
    series
    title
    type
    volume
    year
    eprint
    eprinttype
    eprintclass
    primaryclass
    doi
    url
    lastaccessed
    isbn
    isbn-13
  }
  {}
  { label }

INTEGERS {output.state blank comma colon period}

STRINGS {urlprefix archiveprefix eprinturl doiprefix doiurl isbnprefix isbnurl}

STRINGS { s t t.org u }

STRINGS { p1 p2 p3 page.count }

FUNCTION {dump.stack}
{
    duplicate$ "stack[top] = [" swap$ * "]" * warning$
}

FUNCTION {init.state.consts}
{ #0 'blank  :=
  #1 'comma :=
  #2 'colon :=
  #3 'period :=
}

FUNCTION {init.variables}
{
  "http://arxiv.org/abs/" 'eprinturl := % prefix to make URL from eprint ref
  "arXiv:" 'archiveprefix := % text prefix printed before eprint ref
  "https://doi.org/" 'doiurl := % prefix to make URL from DOI
  "https://isbnsearch.org/isbn/" 'isbnurl := % here you will get more info
  "URL: " 'urlprefix :=
  "DOI: " 'doiprefix :=
  "ISBN: " 'isbnprefix :=
}

% Taken from ACM-Reference-Format.bst

FUNCTION {empty.or.unknown}
{
  %% Examine the top stack entry, and push 1 if it is empty, or
  %% consists only of whitespace, or is a string beginning with two
  %% queries (??), and otherwise, push 0.
  %%
  %% This function provides a replacement for empty$, with the
  %% convenient feature that unknown values marked by two leading
  %% queries are treated the same as missing values, and thus, do not
  %% appear in the output .bbl file, and yet, their presence in .bib
  %% file(s) serves to mark values which are temporarily missing, but
  %% are expected to be filled in eventually once more data is
  %% obtained.  The TeX User Group and BibNet bibliography archives
  %% make extensive use of this practice.
  %%
  %% An empty string cannot serve the same purpose, because just as in
  %% statistics data processing, an unknown value is not the same as an
  %% empty value.
  %%
  %% Before: stack = ... top:[string]
  %% After:  stack = ... top:[0 or 1]

  duplicate$ empty$
    { pop$ #1 }
    { #1 #2 substring$ "??" = }
  if$
}

FUNCTION {writeln}
{
  write$                % output top of the stack
  newline$              % immediate write of newline (not via stack)
}

FUNCTION {blank.mode}
{ blank 'output.state :=
}

FUNCTION {comma.mode}
{ comma 'output.state :=
}

FUNCTION {colon.mode}
{ colon 'output.state :=
}

FUNCTION {period.mode}
{ period 'output.state :=
}

FUNCTION {strip.doi}
{ % Strip any Web address prefix to recover the bare DOI, leaving the
  % result on the output stack, as recommended by CrossRef DOI
  % documentation.
  % For example, reduce "http://doi.acm.org/10.1145/1534530.1534545" to
  % "10.1145/1534530.1534545".  That is later typeset and displayed as
  % DOI: 10.1145/1534530.1534545 as the LAST item in the reference list
  % entry.  Publisher websites wrap this with a suitable link to a real
  % URL to resolve the DOI, and the master https://doi.org/ address is
  % preferred, since publisher-specific URLs can disappear in response
  % to economic events.  All journals are encouraged by the DOI
  % authorities to use that typeset format and link procedures for
  % uniformity across all publications that include DOIs in reference
  % lists.
  % The numeric prefix is guaranteed to start with "10.", so we use
  % that as a test.
  doi #1 #3 substring$ "10." =
    { doi }
    {
      doi 't :=  % get modifiable copy of DOI

      % Change https:// to http:// to strip both prefixes

      t #1 #8 substring$ "https://" =
        { "http://"  t #9 t text.length$ #8 - substring$ * 't := }
        { }
      if$

      t #1 #7 substring$ "http://" =
        {
            t #8 t text.length$ #7 - substring$ 't :=

            "INTERNAL STYLE-FILE ERROR" 's :=

            % search for next "/" and assign its suffix to s

            { t text.length$ }
            {
              t #1 #1 substring$ "/" =
                {
                  % save rest of string as true DOI (should be 10.xxxx/yyyy)
                  t #2 t text.length$ #1 - substring$ 's :=
                  "" 't :=    % empty string t terminates the loop
                }
                {
                  % discard first character and continue loop: t <= substring(t,2,last)
                  t #2 t text.length$ #1 - substring$ 't :=
                }
              if$
            }
            while$

            % check for valid DOI (should be 10.xxxx/yyyy)
            s #1 #3 substring$ "10." =
              { }
              { "unrecognized DOI substring " s * " in DOI value [" * doi * "]" * warning$ }
            if$

            s   % push the stripped DOI on the output stack

        }
        {
          "unrecognized DOI value [" doi * "]" * warning$
          doi   % push the unrecognized original DOI on the output stack
        }
      if$
    }
  if$
}

FUNCTION {output.eprint} 
{ % return with stack empty
  eprint empty.or.unknown
    { }
    { "\showeprint"
         archiveprefix empty.or.unknown
           { eprinttype empty.or.unknown
               { }
               { "[" eprinttype "]" * * * }
             if$
           }
           { "[" archiveprefix "l" change.case$ "]" * * * }
         if$
         "{" *
         primaryclass empty.or.unknown
           { eprintclass empty.or.unknown
             { }
             { eprintclass "/" * * }
             if$
           }
           { primaryclass "/" * * }
         if$
         eprint "}" * *
         writeln
    }
  if$
}

% Do not output url if doi is defined
%
FUNCTION {output.url}
{ % return with stack untouched
  % output URL and associated lastaccessed fields
  doi empty.or.unknown
  {
    url empty.or.unknown
      { }
      {
          %% Use \urldef, outside \showURL, so that %nn, #, etc in URLs work
          %% correctly.  Put the actual URL on its own line to reduce the
          %% likelihood of BibTeX's nasty line wrapping after column 79.
          %% \url{} can undo this, but if that doesn't work for some reason
          %% the .bbl file would have to be repaired manually.
          lastaccessed empty.or.unknown
            { "" }
            { "Retrieved " lastaccessed * " from " * }
          if$
          "\url{" * url * "}" * writeln
      }
      if$
  }
  { }
  if$
}

FUNCTION {output.doi} 
{ % output non-empty DOI as one-line sentence (stack untouched)
  doi empty.or.unknown
    { }
    {
      "\url{" doiurl * strip.doi * "}" * writeln
    }
  if$
}

FUNCTION {output.isbn} 
{ % output non-empty ISBN-10 or ISBN-13 as one-liners (stack untouched)
  %% show both 10- and 13-digit ISBNs
  isbn empty.or.unknown
    { }
    {
      "\url{" isbnurl * isbn * "}" * writeln
    }
  if$
  isbn-13 empty.or.unknown
    { }
    {
      "\url{" isbnurl * isbn-13 * "}" * writeln
    }
  if$
}

FUNCTION {output.nonnull}
{ % stack before: ... r s t  stack after: ... r t   file out: s<punctuation>
  't :=
  output.state comma =
    { ", " * write$ }
    { output.state period =
        { add.period$ writeln }
        { output.state blank =
            { " " * write$ }
            { ": " * write$ }
	  if$
	}
      if$
    }
  if$
  t
}

FUNCTION {output}
{ % discard top token if empty, else like output.nonnull
  duplicate$ empty.or.unknown
    'pop$
    'output.nonnull
  if$
}

FUNCTION {output.check}
{ % like output, but warn if the field on top of the stack is not set
  't :=
  duplicate$ empty.or.unknown
    { pop$ "empty " t * " in " * cite$ * warning$ }
    'output.nonnull
  if$
}

FUNCTION {output.bibitem}
{ newline$
  "\bibitem{" write$
  cite$ write$
  "}" writeln
  ""
  blank 'output.state :=
}

FUNCTION {output.eprint.isbn.doi.url}
{ % enter and return with stack empty
  %% We switch now from buffered output to output of complete lines, so
  %% that the Issue .. URL data have their own lines, and are less likely
  %% to be line-wrapped by BibTeX's short-sighted algorithm, which wraps
  %% lines longer than 79 characters, backtracking to what it thinks is
  %% a break point in the string.  Any such wrapping MUST be undone to
  %% prevent percent-newline from appearing in DOIs and URLs.  The
  %% output data are intentionally wrapped in \showxxx{} macros at
  %% beginning of line, and that supply their own punctuation (if they
  %% are not defined to suppress output entirely), to make it easier for
  %% other software to recover them from .bbl files.
  %%
  %% It also makes it possible to later change the macro definitions
  %% to suppress particular output values, or alter their appearance.

  output.eprint
  output.isbn
  output.doi    % DOI is ALWAYS last according to CrossRef DOI documentation
  output.url    % but ACM wants URL last
}

FUNCTION {output.note}
{ % return with stack empty
  note empty.or.unknown
    { }
    { note add.period$ writeln }
  if$
}

FUNCTION {finish.sentence}
{
  add.period$
  writeln
}

FUNCTION {not}
{   { #0 }
    { #1 }
  if$
}

FUNCTION {and}
{   'skip$
    { pop$ #0 }
  if$
}

FUNCTION {or}
{   { pop$ #1 }
    'skip$
  if$
}

FUNCTION {field.or.null}
{ duplicate$ empty$
    { pop$ "" }
    'skip$
  if$
}

FUNCTION {emphasize}
{ duplicate$ empty$
    { pop$ "" }
    { "\emph{" swap$ * "}" * }
  if$
}

INTEGERS { nameptr namesleft numnames }

FUNCTION {format.names}
{ 's :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
%    { s nameptr "{ff~}{vv~}{ll}{, jj}" format.name$ 't :=
    { s nameptr "{f.~}{vv~}{ll}{, jj}" format.name$ 't :=
%     { s nameptr "{f{.\,}.~}{vv~}{ll}{, jj}" format.name$ 't :=
      nameptr #1 >
        { namesleft #1 >
            { ", " * t * }
            { numnames #2 >
                { "," * }
                'skip$
              if$
              t "others" =
                { " \IfLanguageName{danish}{og andre}{et~al.}" * }
                { " \IfLanguageName{danish}{og}{and} " * t * }
              if$
            }
          if$
        }
        't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}

% For names inside entries (e.g., editors of an "In ...");
% this is exactly abbrv.bst's `format.names' function.
FUNCTION {format.innames}
{ 's :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr "{f.~}{vv~}{ll}{, jj}" format.name$ 't :=
      nameptr #1 >
	{ namesleft #1 >
	    { ", " * t * }
	    { numnames #2 >
		{ "," * }
		'skip$
	      if$
	      t "others" =
		{ " \IfLanguageName{danish}{og andre}{et~al.}" * }
		{ " \IfLanguageName{danish}{og}{and} " * t * }
	      if$
	    }
	  if$
	}
	't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}

FUNCTION {format.authors}
{ author empty.or.unknown
    { "" }
    { author format.names }
  if$
}

FUNCTION {format.editors}
{ editor empty.or.unknown
    { "" }
    { editor format.names 
      editor num.names$ #1 >
        { " (editors)" * }
        { " (editor)" * }
      if$
    }
  if$
}

FUNCTION {format.title}
{ title empty.or.unknown
    { "" }
    { title emphasize }
  if$
}

FUNCTION {lower.title}
{ title empty.or.unknown
    { "" }
    { title "t" change.case$ }
  if$
}

FUNCTION {format.booktitle}
{ booktitle empty.or.unknown
    { "" }
    { booktitle emphasize }
  if$
}

FUNCTION {n.dashify}
{ 't :=
  ""
    { t empty$ not }
    { t #1 #1 substring$ "-" =
	{ t #1 #2 substring$ "--" = not
	    { "--" *
	      t #2 global.max$ substring$ 't :=
	    }
	    {   { t #1 #1 substring$ "-" = }
		{ "-" *
		  t #2 global.max$ substring$ 't :=
		}
	      while$
	    }
	  if$
	}
	{ t #1 #1 substring$ *
	  t #2 global.max$ substring$ 't :=
	}
      if$
    }
  while$
}

FUNCTION {tie.or.space.connect}
{ duplicate$ text.length$ #3 <
    { "~" }
    { " " }
  if$
  swap$ * *
}

FUNCTION {either.or.check}
{ empty$
    'pop$
    { "can't use both " swap$ * " fields in " * cite$ * warning$ }
  if$
}

FUNCTION {format.series}
{ series empty.or.unknown
  { "" }
  { series emphasize
    volume empty$
      'skip$
      { " \textbf{" * volume * "}" * }
    if$
    "volume and number" number either.or.check
  }
  if$
}

FUNCTION {format.edition}
{ edition empty.or.unknown
    { "" }
    { output.state comma =
	{ edition "l" change.case$ "~edition" * }
	{ edition "t" change.case$ "~edition" * }
      if$
    }
  if$
}

FUNCTION {format.pages}
{ pages empty.or.unknown
    { "" }
    { pages n.dashify }
  if$
}

FUNCTION {format.articleno}
{
  articleno empty.or.unknown
    { "" }
    {
      numpages empty.or.unknown
        { "articleno, but no numpages field, in " cite$ * warning$ }
        { }
      if$
    }
  if$
}

FUNCTION {format.page.count}
{
  page.count empty.or.unknown
    { "" }
    { "\bibinfo{numpages}{" page.count * "}~pages" * }
  if$
}

FUNCTION {le}
{
  %% test whether first number is less than or equal to second number
  %% stack before: n1 n2
  %% stack after:  if n1 <= n2 then 1 else 0

  %% "DEBUG: le " cite$ * warning$
  > { #0 } { #1 } if$
}

FUNCTION {is.leading.digit}
{
  %% test whether first character of string is a digit
  %% stack before: string
  %% stack after:  if first-char-is-digit then 1 else 0

  #1 #1 substring$                      % replace string by string[1:1]
  duplicate$                            % string[1:1] string[1:1]
  chr.to.int$
  "0" chr.to.int$ swap$ le              % "0" <= string[1:1] --> 0-or-1
  swap$                                 % 0-or-1 string[1:1]
  chr.to.int$
  "9" chr.to.int$ le                    % string[1:1} <= "9" --> 0-or-1
  and
}

FUNCTION {skip.digits}
{
  duplicate$
  't :=
  't.org :=
  "" 'u :=

  { t text.length$ }
  {
    t is.leading.digit
      { t #2 t text.length$ #1 - substring$ }
      {
        t 'u :=
        ""
      }
    if$
    't :=
  }
  while$

  u                                                      % rest of string
  t.org #1 t.org text.length$ u text.length$ - substring$% leading digits
}

FUNCTION {skip.nondigits}
{
  %% skip over leading nondigits in string
  %% stack before:  string
  %% stack after: rest-of-string

  %% "DEBUG: enter skip.nondigits " cite$ * warning$

  't :=
  "" 'u :=

  { t text.length$ }
  {
    %% "================DEBUG: skip.nondigits   t = [" t * "]" * warning$
    t is.leading.digit
      {
        t 'u :=
        ""
      }
      { t #2 t text.length$ #1 - substring$ }
    if$
    't :=
  }
  while$

  u % rest of string

  %% dump.stack
  %% "DEBUG: leave skip.nondigits " cite$ * warning$
}

FUNCTION {parse.next.number}
{
  %% stack before: string
  %% stack after:  rest-of-string next-numeric-part-of-string
  %% Example:
  %% stack before: "123:1--123:59"
  %% stack after:  ":1--123:59" "123"

  's :=
  s skip.nondigits 's :=
  s skip.digits
}

FUNCTION {reduce.pages.to.page.count}
{
  %% stack before: arbitrary and unused
  %% stack after:  unchanged
  %%
  %% For the new-style pagination with article number and numpages or
  %% pages, we expect to have BibTeX entries containing something like
  %%     articleno = "17",
  %%     pages     = "1--23",
  %% with output "Article 17, 23 pages",
  %% or
  %%     articleno = "17",
  %%     numpages  = "23",
  %% with output "Article 17, 23 pages",
  %% or
  %%     articleno = "17",
  %%     pages     = "17:1--17:23",
  %% with output "Article 17, 23 pages",
  %%
  %% If articleno is missing or empty, then we should output "1--23",
  %% "23" (with a warning of a missing articleno), or "17:1--17:23",
  %% respectively.

  pages
  parse.next.number 'p1 :=
  parse.next.number 'p2 :=
  parse.next.number 'p3 :=
  parse.next.number 'page.count :=

  duplicate$
  empty.or.unknown
    {  }
    {
      duplicate$ "unexpected trailing garbage [" swap$ *
      "] after n:p1--n:p2 in pages = [" *
      pages *
      "] in " *
      cite$ *
      warning$
    }
  if$

  pop$

  p1 p3 =   p2 "1" =   and   numpages empty.or.unknown   and
    { "INFO: reduced pages = [" pages * "] to numpages = [" * page.count * "]" * warning$ }
    {
      numpages empty.or.unknown
        { pages }
        { numpages }
      if$
      'page.count :=
    }
  if$

  p1 "1" =   p3 empty.or.unknown   and   numpages empty.or.unknown   and
    {
      p2 'page.count :=
      "INFO: reduced pages = [" pages * "] to numpages = [" * page.count * "]" * warning$
    }
    {
      numpages empty.or.unknown
        { pages }
        { numpages }
      if$
      'page.count :=
    }
  if$
}

FUNCTION {format.articleno.numpages}
{
  %% There are seven possible outputs, depending on which fields are set.
  %%
  %% These four are handled here:
  %%
  %% articleno, numpages, pages->"Article articleno-value, numpages-value pages"
  %% articleno, numpages -> "Article articleno-value, numpages-value pages"
  %% articleno, pages -> "Article articleno-value, reduced-pages-value pages"
  %% articleno -> "Article articleno-value" and warn about missing numpages
  %%
  %% The remaining three have already been handled by
  %% format.pages.check.without.articleno:
  %%
  %%     numpages, pages -> "pages-value"
  %%     numpages -> "numpages-value"
  %%     pages -> "pages-value"
  %%
  %% We no longer issue warning when missing articleno, but having numpages

  articleno empty.or.unknown 
    { "" }
    {
      numpages empty.or.unknown
        {
          pages empty.or.unknown
            {
              "articleno, but no pages or numpages field in "
                 cite$ * warning$
              "" 'page.count :=
            }
            { reduce.pages.to.page.count }
          if$
        }
        { numpages 'page.count := }
      if$

      %% The Article number is now handled in format.day.month.year because
      %% ACM prefers the style "Digital Libraries 12, 3, Article 5 (July 2008)"
      %% over "Digital Libraries 12, 3 (July 2008), Article 5"
      %% format.articleno output
      format.page.count
    }
  if$
}

FUNCTION {format.pages.check.without.articleno}
{ %% format pages field only if articleno is absent
  %% stack after: pages-specification
  numpages missing$ pages missing$ and
    { "page numbers missing in both pages and numpages fields in " cite$ * warning$ }
    { }
  if$

  articleno empty.or.unknown
    {
      pages missing$
        { numpages }
        { format.pages }
      if$
    }
    { "" }
  if$
}

INTEGERS { multiresult }

FUNCTION {multi.page.check}
{ 't :=
  #0 'multiresult :=
    { multiresult not
      t empty$ not
      and
    }
    { t #1 #1 substring$
      duplicate$ "-" =
      swap$ duplicate$ "," =
      swap$ "+" =
      or or
	{ #1 'multiresult := }
	{ t #2 global.max$ substring$ 't := }
      if$
    }
  while$
  multiresult
}

FUNCTION {format.pp.pages}
{ pages empty$
    { "" }
    { pages multi.page.check
        { "pp.~" pages n.dashify * }
        { "p.~" pages * }
      if$
    }
  if$
}

FUNCTION {format.chapter.pages}
{ chapter empty$
    'format.pages
    { type empty$
        { "Chapter~" chapter * }
        { type "t" change.case$ chapter tie.or.space.connect }
      if$
      pages empty$
        'skip$
        { ", " * format.pages * }
      if$
    }
  if$
}

FUNCTION {empty.misc.check}
{ author empty$ title empty$ howpublished empty$
  month empty$ year empty$ note empty$
  and and and and and
  key empty$ not and
    { "all relevant fields are empty in " cite$ * warning$ }
    'skip$
  if$
}

FUNCTION {format.thesis.type}
{ type empty$
    'skip$
    { pop$
      type
    }
  if$
}

FUNCTION {format.tr.number}
{ type empty$
    { "Technical report" }
    'type
  if$
  number empty$
    { "t" change.case$ }
    { " \textbf{" * number * "}" * }
  if$
}

FUNCTION {format.article.crossref}
{ key empty$
    { journal empty$
	{ "need key or journal for " cite$ * " to crossref " * crossref *
	  warning$
	  ""
	}
	{ "in \emph{" journal * "}" * }
      if$
    }
    { "in " key * }
  if$
  " \cite{" * crossref * "}" *
}

FUNCTION {format.crossref.editor}
{ editor #1 "{vv~}{ll}" format.name$
  editor num.names$ duplicate$
  #2 >
    { pop$ " \IfLanguageName{danish}{og andre}{et~al.}" * }
    { #2 <
	'skip$
	{ editor #2 "{ff }{vv }{ll}{ jj}" format.name$ "others" =
	    { " \IfLanguageName{danish}{og andre}{et~al.}" * }
	    { " \IfLanguageName{danish}{og}{and} " * editor #2 "{vv~}{ll}" format.name$ * }
	  if$
	}
      if$
    }
  if$
}

FUNCTION {format.book.crossref}
{ volume empty$
    { "empty volume in " cite$ * "'s crossref of " * crossref * warning$
      "in "
    }
    { " Vol.~\textbf{" volume * "}" *
      " of " *
    }
  if$
  editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
	{ series empty$
	    { "need editor, key, or series for " cite$ * " to crossref " *
	      crossref * warning$
	      "" *
	    }
%	    { " \emph{" * series * "}" * }
	    { " \emph{" * series * "}" * }
	  if$
	}
	{ key * }
      if$
    }
    { format.crossref.editor * }
  if$
  " \cite{" * crossref * "}" *
}

FUNCTION {format.incoll.inproc.crossref}
{ editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
	{ booktitle empty$
	    { "need editor, key, or booktitle for " cite$ * " to crossref " *
	      crossref * warning$
	      ""
	    }
	    { "in \emph{" booktitle * "}" * }
	  if$
	}
	{ "in " key * }
      if$
    }
    { "in " format.crossref.editor * }
  if$
  " \cite{" * crossref * "}" *
}  

FUNCTION {format.year}
{ year empty.or.unknown
    { month empty$
	    'skip$
	    { "there's a month but no year in " cite$ * warning$ }
      if$
    ""
    }
    { " (" year * ")" * }
  if$
}

FUNCTION {format.isbn}
{ year empty.or.unknown
    { month empty$
	    'skip$
	    { "there's a month but no year in " cite$ * warning$ }
      if$
    ""
    }
    { " (" year * ")" * }
  if$
}

FUNCTION {format.lastaccessed}
{ lastaccessed empty.or.unknown
    { "last-access date not provided in " cite$ * warning$ }
    { lastaccessed * }
  if$
}

FUNCTION {article}
{ output.bibitem
  format.authors "author" output.check
  comma.mode
  lower.title "title" output.check
  crossref missing$
    { journal empty$
        { "empty journal in " cite$ * warning$
          ""
        }
        { journal emphasize
          volume empty$
	    'skip$
	    { " \textbf{" * volume * "}" * }
          if$
          number empty$
	    'skip$
            { ", " * number * }
          if$
        }
      if$
      output
      blank.mode
      format.year "year" output.check
      comma.mode
    }
    { format.article.crossref output.nonnull }
  if$
  format.pages.check.without.articleno output
  format.articleno.numpages output
  finish.sentence
  output.eprint.isbn.doi.url
  output.note
}

FUNCTION {book}
{ output.bibitem
  author empty$
    { editor empty$
        { organization empty$
            { "author, editor, or organization required in " cite$ * warning$ }
            { organization output }
          if$
        }
        { format.editors output }
      if$
    }
    { format.authors output.nonnull
      crossref missing$
	    { "author and editor" editor either.or.check }
	    'skip$
      if$
    }
  if$
  comma.mode
  format.title "title" output.check
  format.series output
  format.edition output
  crossref missing$
    { publisher "publisher" output.check
      address output
    }
    { format.book.crossref output.nonnull
    }
  if$
  blank.mode
  format.year "year" output.check
  finish.sentence
  output.eprint.isbn.doi.url
  output.note
}

FUNCTION {booklet}
{ output.bibitem
  format.authors output
  comma.mode
  format.title "title" output.check
  howpublished output
  address output
  blank.mode
  format.year output
  finish.sentence
  output.eprint.isbn.doi.url
  output.note
}

FUNCTION {inbook}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check }
    { format.authors output.nonnull
      crossref missing$
        { "author and editor" editor either.or.check }
        'skip$
      if$
    }
  if$
  comma.mode
  format.title "title" output.check
  author empty$
    'skip$
    { format.editors output }
  if$
  format.series output
  format.edition output
  crossref missing$
    { publisher "publisher" output.check
      address output
    }
    { format.book.crossref output.nonnull }
  if$
  blank.mode
  format.year "year" output.check
  format.chapter.pages "chapter and pages" output.check
  finish.sentence
  output.eprint.isbn.doi.url
  output.note
}

FUNCTION {incollection}
{ output.bibitem
  format.authors "author" output.check
  comma.mode
  lower.title "title" output.check
  crossref missing$
    { format.booktitle "booktitle" output.check
%      format.editors output
      format.series output      
      format.edition output
      publisher "publisher" output.check
      address output
      blank.mode
      format.year "year" output.check
    }
    { format.incoll.inproc.crossref output.nonnull }
  if$
  comma.mode
  format.chapter.pages output
  finish.sentence
  output.eprint.isbn.doi.url
  output.note
}

FUNCTION {inproceedings}
{ output.bibitem
  format.authors "author" output.check
  comma.mode
  lower.title "title" output.check
  crossref missing$
    { format.booktitle "booktitle" output.check
%      format.editors output
      format.series output
      publisher output
      address output
      blank.mode 
      format.year "year" output.check
    }
    { format.incoll.inproc.crossref output.nonnull }
  if$
  comma.mode
  format.pages output
  finish.sentence
  output.eprint.isbn.doi.url
  output.note
}

FUNCTION {conference} { inproceedings }

FUNCTION {manual}
{ output.bibitem
  author empty$
    { organization output }
    { format.authors output.nonnull }
  if$
  comma.mode
  format.title "title" output.check
  format.edition output
  author empty$
    'skip$
    { organization output }
  if$
  address output
  blank.mode
  format.year output
  finish.sentence
  output.eprint.isbn.doi.url
  output.note
}

FUNCTION {mastersthesis}
{ output.bibitem
  format.authors "author" output.check
  comma.mode
  lower.title "title" output.check
  "M.\,Sc.\ thesis" format.thesis.type output.nonnull
  school "school" output.check
  address output
  blank.mode
  format.year "year" output.check
  finish.sentence
  output.eprint.isbn.doi.url
  output.note
}

FUNCTION {misc}
{ output.bibitem
  format.authors output
  comma.mode
  format.title output
  howpublished output
  blank.mode
  format.year output
  finish.sentence
  output.eprint.isbn.doi.url
  output.note
  empty.misc.check
}

FUNCTION {phdthesis}
{ output.bibitem
  format.authors "author" output.check
  comma.mode
  lower.title "title" output.check
  "Ph.\,D.\ thesis" format.thesis.type output.nonnull
  school "school" output.check
  address output
  blank.mode
  format.year "year" output.check
  finish.sentence
  output.eprint.isbn.doi.url
  output.note
}

FUNCTION {proceedings}
{ output.bibitem
  editor empty$
    { organization output }
    { format.editors output.nonnull }
  if$
  comma.mode
  format.title "title" output.check
  format.editors output
  format.series output
  publisher empty$
    { organization output }
    { publisher output }
  if$
  address output
  blank.mode
  format.year "year" output.check
  finish.sentence
  output.eprint.isbn.doi.url
  output.note
}

FUNCTION {techreport}
{ output.bibitem
  author empty$
    { format.editors "editor" output.check }
    { format.authors output.nonnull }
  if$
  comma.mode
  lower.title "title" output.check
  format.tr.number output.nonnull
  institution "institution" output.check
  address output
  blank.mode
  format.year "year" output.check
  finish.sentence
  output.eprint.isbn.doi.url
  output.note
}

FUNCTION {unpublished}
{ output.bibitem
  format.authors "author" output.check
  comma.mode
  lower.title "title" output.check
  blank.mode
  format.year output
  finish.sentence
  output.eprint.isbn.doi.url
  output.note
}

FUNCTION {website}
{ output.bibitem
  author empty$
    { editor empty$
        { organization empty$
            { "author, editor, or organization required in " cite$ * warning$ }
            { organization output }
          if$
        }
        { format.editors output }
      if$
    }
    { format.authors output }
  if$  
  comma.mode
  lower.title "title" output.check
  "Website " output
  blank.mode
  format.year "year" output.check
  finish.sentence
  output.eprint.isbn.doi.url
  output.note
}

FUNCTION {webpage}
{ output.bibitem
  author empty$
    { editor empty$
        { organization empty$
            { "author, editor, or organization required in " cite$ * warning$ }
            { organization output }
          if$
        }
        { format.editors output }
      if$
    }
    { format.authors output }
  if$  
  comma.mode
  lower.title "title" output.check
  blank.mode
  format.year "year" output.check
  finish.sentence
  output.eprint.isbn.doi.url
  output.note
}

FUNCTION {default.type} { misc }

MACRO {jan} {"Jan."}

MACRO {feb} {"Feb."}

MACRO {mar} {"Mar."}

MACRO {apr} {"Apr."}

MACRO {may} {"May"}

MACRO {jun} {"June"}

MACRO {jul} {"July"}

MACRO {aug} {"Aug."}

MACRO {sep} {"Sept."}

MACRO {oct} {"Oct."}

MACRO {nov} {"Nov."}

MACRO {dec} {"Dec."}

MACRO {acmcs} {"ACM Comput. Surv."}

MACRO {acta} {"Acta Inf."}

MACRO {cacm} {"Commun. ACM"}

MACRO {ibmjrd} {"IBM J. Res. Dev."}

MACRO {ibmsj} {"IBM Syst.~J."}

MACRO {ieeese} {"IEEE Trans. Softw. Eng."}

MACRO {ieeetc} {"IEEE Trans. Comput."}

MACRO {ieeetcad}
 {"IEEE Trans. Comput.-Aided Design Integrated Circuits"}

MACRO {ipl} {"Inf. Process. Lett."}

MACRO {jacm} {"J.~ACM"}

MACRO {jcss} {"J.~Comput. Syst. Sci."}

MACRO {scp} {"Sci. Comput. Programming"}

MACRO {sicomp} {"SIAM J. Comput."}

MACRO {tocs} {"ACM Trans. Comput. Syst."}

MACRO {tods} {"ACM Trans. Database Syst."}

MACRO {tog} {"ACM Trans. Gr."}

MACRO {toms} {"ACM Trans. Math. Softw."}

MACRO {toois} {"ACM Trans. Office Inf. Syst."}

MACRO {toplas} {"ACM Trans. Program. Lang. Syst."}

MACRO {tcs} {"Theoretical Comput. Sci."}

READ

FUNCTION {sortify}
{ purify$
  "l" change.case$
}

INTEGERS { len }

FUNCTION {chop.word}
{ 's :=
  'len :=
  s #1 len substring$ =
    { s len #1 + global.max$ substring$ }
    's
  if$
}

FUNCTION {sort.format.names}
{ 's :=
  #1 'nameptr :=
  ""
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { nameptr #1 >
	{ "   " * }
	'skip$
      if$
      s nameptr "{ll{ }}{  f{ }}{  jj{ }}" format.name$ 't :=
      nameptr numnames = t "others" = and
	{ "\IfLanguageName{danish}{og andre}{et~al}" * }
	{ t sortify * }
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}

FUNCTION {sort.format.title}
{ 't :=
  "A " #2
    "An " #3
      "The " #4 t chop.word
    chop.word
  chop.word
  sortify
  #1 global.max$ substring$
}

FUNCTION {author.sort}
{ author empty$
    { key empty$
	{ "to sort, need author or key in " cite$ * warning$
	  ""
	}
	{ key sortify }
      if$
    }
    { author sort.format.names }
  if$
}

FUNCTION {author.editor.sort}
{ author empty$
    { editor empty$
	{ key empty$
	    { "to sort, need author, editor, or key in " cite$ * warning$
	      ""
	    }
	    { key sortify }
	  if$
	}
	{ editor sort.format.names }
      if$
    }
    { author sort.format.names }
  if$
}

FUNCTION {author.organization.sort}
{ author empty$
    { organization empty$
	{ key empty$
	    { "to sort, need author, organization, or key in " cite$ * warning$
	      ""
	    }
	    { key sortify }
	  if$
	}
	{ "The " #4 organization chop.word sortify }
      if$
    }
    { author sort.format.names }
  if$
}

FUNCTION {editor.organization.sort}
{ editor empty$
    { organization empty$
	{ key empty$
	    { "to sort, need editor, organization, or key in " cite$ * warning$
	      ""
	    }
	    { key sortify }
	  if$
	}
	{ "The " #4 organization chop.word sortify }
      if$
    }
    { editor sort.format.names }
  if$
}

FUNCTION {presort}
{ type$ "book" =
  type$ "inbook" =
  or
    'author.editor.sort
    { type$ "proceedings" =
	'editor.organization.sort
	{ type$ "manual" =
	    'author.organization.sort
	    'author.editor.sort
	  if$
	}
      if$
    }
  if$
  "    "
  *
  year field.or.null sortify
  *
  "    "
  *
  title field.or.null
  sort.format.title
  *
  #1 entry.max$ substring$
  'sort.key$ :=
}

ITERATE {presort}

SORT

STRINGS { longest.label }

INTEGERS { number.label longest.label.width }

FUNCTION {initialize.longest.label}
{ "" 'longest.label :=
  #1 'number.label :=
  #0 'longest.label.width :=
}

FUNCTION {longest.label.pass}
{ number.label int.to.str$ 'label :=
  number.label #1 + 'number.label :=
  label width$ longest.label.width >
    { label 'longest.label :=
      label width$ 'longest.label.width :=
    }
    'skip$
  if$
}

EXECUTE {initialize.longest.label}

ITERATE {longest.label.pass}

FUNCTION {begin.bib}
{ preamble$ empty$
    'skip$
    { preamble$ write$ newline$ }
  if$
  "\begin{thebibliography}{"  longest.label  * "}" * writeln
  "\ifx \showDOI      \undefined \def \showDOI       #1{#1}\fi" writeln
  "\ifx \showISBNx    \undefined \def \showISBNx     #1{#1}\fi" writeln
  "\ifx \showISBNxiii \undefined \def \showISBNxiii  #1{#1}\fi" writeln
  "\ifx \showURL      \undefined \def \showURL       #1{#1}\fi" writeln
  "% The following commands are used for tagged output and should be " writeln
  "% invisible to TeX" writeln
  "\providecommand\bibinfo[2]{#2}" writeln
  "\providecommand\showeprint[2][]{arXiv:#2}" writeln
}

EXECUTE {begin.bib}
EXECUTE {init.variables}
EXECUTE {init.state.consts}

ITERATE {call.type$}

FUNCTION {end.bib}
{ newline$
  "\end{thebibliography}" writeln
}

EXECUTE {end.bib}


