%%
%% This is the file `degruyter-journal-a.cls'.
%%
\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{degruyter-journal-a}[2012/05/22 v1.4]
%%
%% LaTeX document class for De Gruyter journals, format 155x230.
%% By Dimler & Albroscheit Partnerschaft (www.dimler-albroscheit.de).
%%
%% v1.0 -- Welcome!
%% v1.1 -- J.~Group~Theory; J.~Inverse~Ill-Posed~Probl; new spacing before \paragraph;
%%         new scaling factor for helvet; \dedicated replaces \dedication
%% v1.2 -- Adv.~Nonlinear~Anal.; keywords and classification ragged right;
%%         theorem environment `notation'
%% v1.3 -- New journal information including DOI abbreviations;
%%         packages `caption' and `subfig'; theorem environment `problem';
%%         new title length for `jnetdy'
%% v1.4 -- Class options `conference' and `bookreview'
%%
\RequirePackage{ifthen}
%
\newboolean{acv}        % Adv. Calc. Var.
\newboolean{anona}      % Adv. Nonlinear Anal.
\newboolean{apam}       % Adv. Pure Appl. Math.
\newboolean{eqc}        % Economic Quality Control
\newboolean{forum}      % Forum Math.
\newboolean{gcc}        % Groups Complex. Cryptol.
\newboolean{gmj}        % Georgian Math. J.
\newboolean{integers}   % Integers
\newboolean{jaa}        % Journal of Applied Analysis
\newboolean{jgt}        % J. Group Theory
\newboolean{jisys}      % J. Intell. Syst.
\newboolean{jip}        % J. Inverse Ill-Posed Probl.
\newboolean{jmc}        % J. Math. Cryptol.
\newboolean{jnetdy}     % J. Non-Equilib. Thermodyn.
\newboolean{mcma}       % Monte Carlo Methods Appl.
\newboolean{rose}       % Random Oper. Stoch. Equ.
%
\newboolean{ak}         % Galley proofs
\newboolean{aop}        % ahead-of-print
\newboolean{communicatedby} % Naming the communicating editor
\newboolean{seceqcount} % equations numbered by section
\newboolean{games}      % games section (INTEGERS)
\newboolean{MTPro}      % mtpro2 fonts
\newboolean{nothanks}   % no researchsupported (use \makethankshere instead)
\newboolean{online}     % --- no functionality ---
\newboolean{noheadline} % no headline, no author
\newboolean{prefaceetc} % for prefaces etc.
\newboolean{noauthor}   % for articles without author (to be used together with "prefaceetc")
\newboolean{nohistory}  % for articles without article history (to be used together with "prefaceetc")
\newboolean{conference} % for conference announcements (calls prefaceetc,noauthor,nohistory)

\newif\if@restonecol
\newif\if@titlepage
\@titlepagefalse

\DeclareOption{landscape}{%
	\setlength{\@tempdima}{\paperheight}%
	\setlength{\paperheight}{\paperwidth}%
	\setlength{\paperwidth}{\@tempdima}
}

\DeclareOption{acv}{\setboolean{acv}{true}}
\DeclareOption{anona}{\setboolean{anona}{true}}
\DeclareOption{apam}{\setboolean{apam}{true}}
\DeclareOption{forum}{\setboolean{forum}{true}}
\DeclareOption{gmj}{\setboolean{gmj}{true}}
\DeclareOption{integers}{\setboolean{integers}{true}}
\DeclareOption{jgt}{\setboolean{jgt}{true}}
\DeclareOption{jaa}{\setboolean{jaa}{true}}
\DeclareOption{jip}{\setboolean{jip}{true}}
\DeclareOption{jmc}{\setboolean{jmc}{true}}
\DeclareOption{jnetdy}{\setboolean{jnetdy}{true}}
\DeclareOption{mcma}{\setboolean{mcma}{true}}
\DeclareOption{rose}{\setboolean{rose}{true}}
\DeclareOption{gcc}{\setboolean{gcc}{true}}
\DeclareOption{eqc}{\setboolean{eqc}{true}}
\DeclareOption{jisys}{\setboolean{jisys}{true}}

\DeclareOption{ak}{\setboolean{ak}{true}}
\DeclareOption{aop}{\setboolean{aop}{true}}
\DeclareOption{games}{\setboolean{games}{true}}
\DeclareOption{MTPro}{\setboolean{MTPro}{true}}
\DeclareOption{nothanks}{\setboolean{nothanks}{true}}
\DeclareOption{online}{\setboolean{online}{true}}
\DeclareOption{noheadline}{\setboolean{noheadline}{true}}
\DeclareOption{prefaceetc}{\setboolean{prefaceetc}{true}}
\DeclareOption{noauthor}{\setboolean{noauthor}{true}}
\DeclareOption{nohistory}{\setboolean{nohistory}{true}}
\DeclareOption{conference}{%
  \setboolean{prefaceetc}{true}%
  \setboolean{noauthor}{true}%
  \setboolean{nohistory}{true}}
\DeclareOption{bookreview}{} %% Now relevance in LaTeX. Just an indicator.
%
\DeclareOption{oneside}{\@twosidefalse \@mparswitchfalse}

\DeclareOption{twoside}{\@twosidetrue  \@mparswitchtrue}
\DeclareOption{draft}{\setlength\overfullrule{5pt}}

\DeclareOption{final}{\setlength\overfullrule{0pt}}
\DeclareOption{titlepage}{\@titlepagetrue}
\DeclareOption{notitlepage}{\@titlepagefalse}
\DeclareOption{onecolumn}{\@twocolumnfalse}

\DeclareOption{twocolumn}{\@twocolumntrue}
\DeclareOption{leqno}{\input{leqno.clo}}
\DeclareOption{fleqn}{\input{fleqn.clo}}

\ifthenelse{\boolean{jnetdy}}%
	{\ExecuteOptions{onecolumn,12pt,final}}%
	{\ExecuteOptions{onecolumn,11pt,final}}%

\ProcessOptions

\ifthenelse{\boolean{jnetdy}}%
	{\input{bk12.clo}}%
	{\input{bk11.clo}}%

%---------------------%
% Changes to bk11.clo %
%---------------------%

\def\itemize{%
	\ifnum \@itemdepth >\thr@@\@toodeep%
	\else
		\advance\@itemdepth\@ne
		\edef\@itemitem{labelitem\romannumeral\the\@itemdepth}%
		\expandafter
		\list
			\csname\@itemitem\endcsname
			{\def\makelabel##1{\hss\llap{\normalfont ##1}}}% % make roman labels upright
	\fi}

\def\@listi{\leftmargin\leftmargini
	\parsep 2.5\p@ \@plus2\p@ \@minus\p@  %NEW (instead of 4.5 + 2.0 - 1.0)
	\topsep 6\p@ \@plus2\p@ \@minus2\p@   %NEW (instead of 9.0 + 3.0 - 5.0)
	\itemsep 3\p@ \@plus2\p@ \@minus\p@   %NEW (instead of 4.5 + 2.0 - 1.0)
}

\let\@listI\@listi
\@listi

\def\@listii{%
	\leftmargin\leftmarginii
	\labelwidth\leftmarginii
	\advance\labelwidth-\labelsep
	\topsep    2.5\p@ \@plus2\p@ \@minus\p@  %NEW (instead of 4.5 + 2.0 - 1.0)
	\parsep    2\p@ \@plus\p@  \@minus\p@
	\itemsep   \parsep
}
%---------------------%


%---------------------%
% make refs allways upright
%---------------------%
\def\@setref#1#2#3{%
	\ifx#1\relax
		\protect\G@refundefinedtrue
		\nfss@text{\reset@font\bfseries ??}%
		\@latex@warning{Reference `#3' on page \thepage \space undefined}%
	\else
		{\textup{\expandafter#2#1\null}}%        % normalfont added
	\fi
}

\def\ref#1{%
	\expandafter\@setref\csname r@#1\endcsname\@firstoftwo{#1}%
}

%---------------------%



%--------------------------------------------------------------------------%
% TYPE AREA                                                                %
% Use class option `dvips' if you compile with latex rather than pdflatex. %
%--------------------------------------------------------------------------%

\usepackage[%
  mag=954,%     % de Gruyter magic number
  asymmetric,%  % margins are not swapped on alternate pages
  hmarginratio=1:1,%
  paperwidth=155truemm,
  paperheight=230truemm,
  textwidth=120truemm,
  textheight=180truemm,
  lines=40,
  top=22truemm,
  headheight=9.5truemm,
  headsep=5.7truemm,
  bindingoffset=0mm]{geometry}[2008/12/21]
%% Why do we need version 4.2 of the geometry package?
%% --> Fixed pdfTEX magnification problem.

\marginparwidth 11truemm%

%--------------------------------------------------------------------------%

\parindent 4.2mm

\renewcommand{\baselinestretch}{1.0}
\renewcommand{\arraystretch}{1.2}
\setcounter{secnumdepth}{2}

\setlength{\jot}{5pt}
\setlength{\lineskip}{1\p@}
\setlength{\normallineskip}{1\p@}
\setlength{\parskip}{0\p@ \@plus \p@}
\setlength{\columnsep}{10\p@}
\setlength{\columnseprule}{0\p@}

\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301

\setcounter{topnumber}{5}              % Max. Anzahl Gleitobjekte oben
\renewcommand{\topfraction}{0.8}       % Max. Anteil Gleitobjekte oben
\setcounter{bottomnumber}{5}           % Max. Anzahl Gleitobjekte unten
\renewcommand{\bottomfraction}{0.8}    % Max. Anteil Gleitobjekte unten
\setcounter{totalnumber}{10}           % Max. Anzahl Gleitobjekte pro Kolumne
\renewcommand{\textfraction}{0.2}      % Min. Anteil Fließtext pro Kolumne
\renewcommand{\floatpagefraction}{0.2} % Min. Anteil Gleitobjekte pro Gleitobjektkolumne (die "[p]-Seiten")
\setcounter{dbltopnumber}{2}
\renewcommand{\dbltopfraction}{0.7}
\renewcommand{\dblfloatpagefraction}{0.5}

\setlength{\@fptop}{0pt}    % Stretch at top of float page
                            % (Lässt p-Objekte oben stehen.)
\setlength{\@dblfptop}{0pt} % Analog für "double-column figures"

\setlength{\floatsep}{16pt plus 2pt minus 2pt} % Abstand zwischen t-Objekten etwas
                                                % größer als Standard (12+2-2)
\setlength{\dblfloatsep}{\floatsep}              % Analog für "double-column figures"

\setlength{\@fpsep}{\floatsep}       % Space between floats on float page (Abstand zwischen
                                     % p-Objekten entspricht Abstand zwischen t-Objekten.)
\setlength{\@dblfpsep}{\dblfloatsep} % Analog für "double-column figures"

%---------%
%\mag=954 %
%---------%

\usepackage{amsmath}
\usepackage{amsthm}[2004/08/06] %Find changes to this package below.
\usepackage{calc}
\usepackage{color}
\usepackage{delarray}
\usepackage{enumerate} %Find changes to this package below.
\usepackage{epsfig}
\usepackage{graphicx}
\usepackage{latexsym}
\usepackage{multicol}
\usepackage{pifont}
\usepackage{rawfonts}
\usepackage[figuresright]{rotating}
\usepackage{url}
\usepackage{verbatim}
\usepackage{xspace}
\usepackage{marginnote}
\usepackage[noadjust]{cite}
\usepackage{ragged2e}

\def\@cite#1#2{\leavevmode \cite@adjust{\upshape
  \citeleft{#1\if@tempswa \citemid #2\fi
  \spacefactor\@m % punctuation in note doesn't affect outside
  }\citeright}}

%------------%
% Hyperlinks %
%------------%
% wenn implicit=true (default), muß \maketitle und \appendix, mit \renewcommand definiert werden.
%\usepackage[implicit=false,bookmarks=false,colorlinks,urlcolor=black,a4paper=false,setpagesize=false]{hyperref}
\usepackage[%
  implicit=false,%
  bookmarks=false,
  pdfpagelabels=false,%
  colorlinks,%
  urlcolor=black,%
  a4paper=false,%
  setpagesize=false%
]{hyperref}
%------------%

%-------%
% FONTS %
%-------%
\usepackage[T1]{fontenc}
\usepackage{textcomp}
\usepackage[scaled=0.88]{helvet} % = 0.92 * \mag
%\usepackage[mathscr]{eucal}

\renewcommand{\rmdefault}{ptm}  %= Adobe Times
\renewcommand{\ttdefault}{cmtt} %= Computer Modern Typewriter

\ifthenelse{\boolean{MTPro}}%
	{\usepackage[mtphrb,mtpfrak,mtpcal,straightbraces,subscriptcorrection]{mtpro2}[2009/1/30]}
	{\usepackage{amsfonts,amssymb}}
%% The packages AMSFONTS, AMSSYMB, EXSCALE are not compatible with mtpro2!!!
%-------%

\newcounter{numberauthors}
%
\newcommand\authorone[1]{\def\@authorone{#1}}
\newcommand\lastnameone[1]{\def\@lastnameone{#1}}
\newcommand\firstnameone[1]{\def\@firstnameone{#1}}
\newcommand\nameshortone[1]{\def\@nameshortone{#1}}
\newcommand\addressone[1]{\def\@addressone{#1}}
\newcommand\countryone[1]{\def\@countryone{#1}}
\newcommand\emailone[1]{\def\@emailone{#1}}
%
\newcommand\authortwo[1]{\def\@authortwo{#1}}
\newcommand\lastnametwo[1]{\def\@lastnametwo{#1}}
\newcommand\firstnametwo[1]{\def\@firstnametwo{#1}}
\newcommand\nameshorttwo[1]{\def\@nameshorttwo{#1}}
\newcommand\addresstwo[1]{\def\@addresstwo{#1}}
\newcommand\countrytwo[1]{\def\@countrytwo{#1}}
\newcommand\emailtwo[1]{\def\@emailtwo{#1}}
%
\newcommand\authorthree[1]{\def\@authorthree{#1}}
\newcommand\lastnamethree[1]{\def\@lastnamethree{#1}}
\newcommand\firstnamethree[1]{\def\@firstnamethree{#1}}
\newcommand\nameshortthree[1]{\def\@nameshortthree{#1}}
\newcommand\addressthree[1]{\def\@addressthree{#1}}
\newcommand\countrythree[1]{\def\@countrythree{#1}}
\newcommand\emailthree[1]{\def\@emailthree{#1}}
%
\newcommand\authorfour[1]{\def\@authorfour{#1}}
\newcommand\lastnamefour[1]{\def\@lastnamefour{#1}}
\newcommand\firstnamefour[1]{\def\@firstnamefour{#1}}
\newcommand\nameshortfour[1]{\def\@nameshortfour{#1}}
\newcommand\addressfour[1]{\def\@addressfour{#1}}
\newcommand\countryfour[1]{\def\@countryfour{#1}}
\newcommand\emailfour[1]{\def\@emailfour{#1}}
%
\newcommand\authorfive[1]{\def\@authorfive{#1}}
\newcommand\lastnamefive[1]{\def\@lastnamefive{#1}}
\newcommand\firstnamefive[1]{\def\@firstnamefive{#1}}
\newcommand\nameshortfive[1]{\def\@nameshortfive{#1}}
\newcommand\addressfive[1]{\def\@addressfive{#1}}
\newcommand\countryfive[1]{\def\@countryfive{#1}}
\newcommand\emailfive[1]{\def\@emailfive{#1}}
%
\newcommand\headlinetitle[1]{\def\@headlinetitle{#1}}
\newcommand\headlineauthor[1]{\def\@shortheadlineauthor{#1}}
\newcommand\abstract[1]{\def\@abstract{#1}}
\newcommand\keywords[1]{\def\@keywords{#1}}
\newcommand\classificationscheme[1]{\def\@classificationscheme{#1}}
\newcommand\classification[1]{\def\@classification{#1}}
\newcommand\communicated[1]{\def\@communicated{#1}}
\newcommand\dedicated[1]{\def\@dedicated{#1}}
\newcommand\acknowledgments[1]{\def\@acknowledgments{#1}}
\newcommand\researchsupported[1]{\def\@researchsupported{#1}}
\newcommand\received[1]{\def\@received{#1}}
\newcommand\revised[1]{\def\@revised{#1}}
\newcommand\accepted[1]{\def\@accepted{#1}}
\newcommand\publishedonline[1]{\def\@publishedonline{#1}}
\newcommand\pubyear[1]{\def\@pubyear{#1}}
\newcommand\aopyear[1]{\def\@aopyear{#1}}
\newcommand\journalname[1]{\def\@journalname{#1}}
\newcommand\journalnameshort[1]{\def\@journalnameshort{#1}}
\newcommand\volume[1]{\def\@volume{#1}}
\newcommand\doipref[1]{\def\@doipref{#1}}
\newcommand\doisep[1]{\def\@doisep{#1}}
\newcommand\doi[1]{\def\@doi{#1}}
\newcommand\doiyear[1]{\def\@doiyear{#1}}

\newcounter{startpage}

\ifthenelse{\boolean{aop}}%
	{\newcommand\firstpage[1]{\setcounter{startpage}{1}}}
	{%
		\ifthenelse{\boolean{ak}}%
			{\newcommand\firstpage[1]{\setcounter{startpage}{1}}}
			{\newcommand\firstpage[1]{\setcounter{startpage}{#1}}}
	}

%----------%
% Defaults %
%----------%

\lastnameone{}
\firstnameone{}
\authorone{\mbox{\@firstnameone~\@lastnameone}}
\nameshortone{\@authorone}
\emailone{}

\lastnametwo{}
\firstnametwo{}
\authortwo{\mbox{\@firstnametwo~\@lastnametwo}}
\nameshorttwo{\@authortwo}
\emailtwo{}

\lastnamethree{}
\firstnamethree{}
\authorthree{\mbox{\@firstnamethree~\@lastnamethree}}
\nameshortthree{\@authorthree}
\emailthree{}

\lastnamefour{}
\firstnamefour{}
\authorfour{\mbox{\@firstnamefour~\@lastnamefour}}
\nameshortfour{\@authorfour}
\emailfour{}

\lastnamefive{}
\firstnamefive{}
\authorfive{\mbox{\@firstnamefive~\@lastnamefive}}
\nameshortfive{\@authorfive}
\emailfive{}

\volume{\textbf{?}}
\pubyear{\textbf{????}}
\aopyear{\@pubyear}
\doiyear{\@pubyear}
\doi{\textbf{???}}
\firstpage{1}
\received{\textbf{???}}
\revised{}
\accepted{}
\publishedonline{}
\headlinetitle{}
\headlineauthor{}
\abstract{}
\keywords{}
\classification{}
\acknowledgments{}
\researchsupported{}
\dedicated{}

%---------------------%
% Journal information %
%---------------------%

%% no journal selected
\journalname{\mbox{}}
\journalnameshort{JOURNAL}
\doipref{JOURNAL}
\doisep{-}
\classificationscheme{MSC}

\newcommand\journaldata[8][\mbox{}]{
  \ifthenelse{\boolean{#2}}{        % Boolean
    \journalname{#1}                % Journal name (forum only)
    \journalnameshort{#3}           % Official abbreviation
    \doipref{#4}                    % DOI abbreviation
    \doisep{#5}                     % Either . (old) or - (new)
    \classificationscheme{#6}       % MSC or PACS
    \setboolean{communicatedby}{#7} % Whether the communicating editor is mentioned
    \setboolean{seceqcount}{#8}}{}} % Whether equations are numbered as (1.1) etc.

\journaldata{acv}{Adv.~Calc.\,Var.}{acv}{-}{MSC}{true}{true}
\journaldata{anona}{Adv.~Nonlinear~Anal.}{anona}{-}{MSC}{false}{true}
\journaldata{apam}{Adv.~Pure~Appl.~Math.}{apam}{-}{MSC}{false}{true}
\journaldata[Forum~Mathematicum]{forum}{Forum~Math.}{forum}{-}{MSC}{true}{true}
\journaldata{gmj}{Georgian~Math.~J.}{gmj}{-}{MSC}{false}{true}
\journaldata{eqc}{Economic~Quality~Control}{eqc}{-}{MSC}{false}{false}
\journaldata{gcc}{Groups~Complex.~Cryptol.}{gcc}{-}{MSC}{false}{false}
\journaldata{integers}{Integers}{integers}{-}{MSC}{false}{false}
\journaldata{jgt}{J.~Group Theory}{jgt}{-}{MSC}{true}{true}
\journaldata{jaa}{Journal~of~Applied~Analysis}{jaa}{-}{MSC}{false}{true}
\journaldata{jip}{J.~Inverse~Ill-Posed~Probl.}{jip}{-}{MSC}{false}{true}
\journaldata{jisys}{J.~Intell.~Syst.}{jisys}{-}{MSC}{false}{false}
\journaldata{jmc}{J.~Math.~Cryptol.}{jmc}{-}{MSC}{true}{true}
\journaldata{jnetdy}{J.~Non-Equilib.~Thermodyn.}{jnetdy}{-}{PACS}{true}{false}
\journaldata{mcma}{Monte~Carlo~Methods~Appl.}{mcma}{-}{MSC}{false}{true}
\journaldata{rose}{Random~Oper.~Stoch.~Equ.}{rose}{-}{MSC}{true}{true}

%---------------%
% More Defaults %
%---------------%

\ifthenelse{\boolean{communicatedby}}
  {\communicated{???}}
  {\communicated{}}

\ifthenelse{\equal{\@classificationscheme}{MSC}}
{\newcommand\classificationname{2010 Mathematics Subject Classification}}
{\newcommand\classificationname{PACS\textsuperscript{\textregistered}(2010)}}

% -------------------------------------------------------------------------------- %
% Header
% -------------------------------------------------------------------------------- %

\usepackage{fancyhdr}
\pagestyle{fancy}
\fancyhead{}\fancyfoot{}
\fancyhead[RO,LE]{\normalsize\thepage}
\fancyhead[CO]{\small\rightmark} %% \rightmark is defined in \def\makeheadlines{...}
\fancyhead[CE]{\small\leftmark}  %% \leftmark is defined in \def\makeheadlines{...}

\renewcommand\sectionmark[1]{}
\renewcommand\subsectionmark[1]{}
\fancypagestyle{mytitlepage}{%
	\fancyhf{}%
	\lhead{%
		\small%
		\ifthenelse{\boolean{aop}}%
			{\@journalnameshort, Ahead of Print\\}
			{%
				\ifthenelse{\boolean{ak}}%
{\@journalnameshort~\textbf{\textcolor{red}{Galley Proof}},~\pageref{LastPage} pages\\}
{\@journalnameshort~\ifthenelse{\boolean{jnetdy}}{\@volume}{\textbf{\@volume}}~(\@pubyear),~\pagerange\\}%
			}%
		DOI 10.1515{\kern1pt}/{\kern1pt}\@doipref\@doisep\@doiyear\@doisep\@doi%
	}%
	
	\ifthenelse{\boolean{games}}
		{\chead{\small\textit{Games Article}}}%
		{}%
	
	\rhead{%
		\small%
		\textbf{\@journalname}\\
		\copyright~de~Gruyter~%
		\ifthenelse{\boolean{aop}}
			{\@aopyear}
			{\@pubyear}
	}
}

%\renewcommand\headrulewidth{30pt}

% -------------------------------------------------------------------------------- %

\newcommand{\pagerange}{%
	\ifthenelse{\equal{\thestartpage}{\pageref{LastPage}}}%
		{\thestartpage}%
		{\thestartpage{\kern1pt}--{\kern1pt}\pageref{LastPage}\ }% %don't delete the last \
}

%----------------------%

\newcommand{\makeauthor}{
	%
	% First page
	%
	\ifthenelse{\equal{\@lastnametwo}{}}%
		{%
			\setcounter{numberauthors}{1}%
			\author{\xprint{\@authorone}}%
		}
		{%
			\ifthenelse{\equal{\@lastnamethree}{}}%
				{%
					\setcounter{numberauthors}{2}%
					\author{\xprint{\@authorone} and \xprint{\@authortwo}}
				}
				{%
					\ifthenelse{\equal{\@lastnamefour}{}}%
						{%
							\setcounter{numberauthors}{3}%
							\author{\xprint{\@authorone}, \xprint{\@authortwo} and \xprint{\@authorthree}}%
						}
						{%
							\ifthenelse{\equal{\@lastnamefive}{}}%
								{%
									\setcounter{numberauthors}{4}%
									\author{\xprint{\@authorone}, \xprint{\@authortwo}, \xprint{\@authorthree} and \xprint{\@authorfour}}%
								}
								{%
									\setcounter{numberauthors}{5}%
									\author{\xprint{\@authorone}, \xprint{\@authortwo}, \xprint{\@authorthree}, \xprint{\@authorfour} and \xprint{\@authorfive}}%
								}
						}%
				}
		}
	%
	% Running title
	%
	\ifthenelse{\equal{\@shortheadlineauthor}{}}%
		{\ifthenelse{\equal{\value{numberauthors}}{1}}%
			{\def\@headlineauthor{\xprint{\@nameshortone}}}%
			{\ifthenelse{\equal{\value{numberauthors}}{2}}%
				{\def\@headlineauthor{\xprint{\@nameshortone} and \xprint{\@nameshorttwo}}}%
				{\ifthenelse{\equal{\value{numberauthors}}{3}}%
					{\def\@headlineauthor{\xprint{\@nameshortone}, \xprint{\@nameshorttwo} and \xprint{\@nameshortthree}}}%
					{\ifthenelse{\equal{\value{numberauthors}}{4}}%
						{\def\@headlineauthor{\xprint{\@nameshortone}, \xprint{\@nameshorttwo}, \xprint{\@nameshortthree} and \xprint{\@nameshortfour}}}
						{\def\@headlineauthor{\xprint{\@nameshortone}, \xprint{\@nameshorttwo}, \xprint{\@nameshortthree}, \xprint{\@nameshortfour} and \xprint{\@nameshortfive}}}%
					}%
				}%
			}%
		}
		{\def\@headlineauthor{\@shortheadlineauthor}}
} % of \newcommand{\makeauthor}

\newdimen\runninglength
\runninglength=100mm

\newdimen\emptylength
\emptylength=1mm

\newcounter{titletoolong}
\setcounter{titletoolong}{0}
\newcounter{authortoolong}
\setcounter{authortoolong}{0}

\def\makeheadlines{%
	%
	% Authors
	%
	\newbox\authrun
	\global\setbox\authrun=\hbox{\small\rmfamily\unboldmath\ignorespaces\@headlineauthor\unskip}%
	\ifdim\wd\authrun<\runninglength
		% okay
	\else
		\setcounter{authortoolong}{1}
		\global\setbox\authrun=\hbox{\small\rmfamily
		\hspace{10mm} Author names too long for running head. Supply a shorter form with \string\headlineauthor}%
		\xdef\@headlineauthor{\copy\authrun}%
	\fi
	%
	% Title
	%
	\newbox\titrun
	\newbox\shorttitrun
	\global\setbox\titrun=\hbox{\small\rmfamily\unboldmath\ignorespaces\@title\unskip}%
	\global\setbox\shorttitrun=\hbox{\small \rmfamily\unboldmath\ignorespaces\@headlinetitle\unskip}%
	%
	\ifdim\wd\titrun>\runninglength
		\setcounter{titletoolong}{1}
		%
		\ifdim\wd\shorttitrun<\emptylength
			\global\setbox\titrun=\hbox{\small\rmfamily Title too long for running head. Supply a shorter form with \string\headlinetitle}
			\xdef\@headlinetitletwo{\copy\titrun}%
			\markboth{\@headlineauthor}{\@headlinetitletwo}%
		\else
			\ifdim\wd\shorttitrun>\runninglength
				\global\setbox\titrun=\hbox{\small\rmfamily Title too long for running head. Supply a shorter form with \string\headlinetitle}%
				\xdef\@headlinetitletwo{\copy\titrun}%
				\markboth{\@headlineauthor}{\@headlinetitletwo}%
			\else
				\setcounter{titletoolong}{0}
				\markboth{\@headlineauthor}{\@headlinetitle}%
			\fi
		\fi
	\else
		 \ifdim\wd\shorttitrun<\emptylength
			\markboth{\@headlineauthor}{\@title}%
		\else
			\markboth{\@headlineauthor}{\@headlinetitle}%
		\fi
	\fi
} % of \def\makeheadlines

\newcommand{\maketitle}{%
	\par
	\begingroup
	\makeauthor%
	\makeheadlines%
	\renewcommand{\thefootnote}{\@fnsymbol\c@footnote}%
	\def\@makefnmark{\rlap{\@textsuperscript{\normalfont\@thefnmark}}}%
	\long\def\@makefntext##1{%
		\parindent 1em\noindent
		\hb@xt@1.8em{\hss\@textsuperscript{\normalfont\@thefnmark}}##1%
	}%
	%
	\if@twocolumn
		\ifnum \col@number=\@ne
			\@maketitle
		\else
			\twocolumn[\@maketitle]%
		\fi
	\else
		\newpage
		\global\@topnum\z@   % Prevents figures from going at top of page.
		\@maketitle
	\fi
	\@thanks
	\endgroup
	\setcounter{footnote}{0}%
	\global\let\thanks\relax
	\global\let\maketitle\relax
	\global\let\@maketitle\relax
	\global\let\@thanks\@empty
	\global\let\@author\@empty
	\global\let\@date\@empty
	\global\let\@title\@empty
	\global\let\title\relax
	\global\let\author\relax
	\global\let\date\relax
	\global\let\and\relax
} % of \newcommand{\maketitle}{%

%-------------------%
% researchsupported %
%-------------------%
\newcommand{\mymakefntext}{%
	\long\def\@makefntext##1{%
	\parindent 0em% %instead of 1em
	\parbox[t]{.04\textwidth}{\hfill\@makefnmark}%
	\hfill%
	\parbox[t]{.95\textwidth}{##1}}%
}

\mymakefntext

\newcommand{\makethanks}[1]{% Don't delete the %-sign at the end of each line!
	\let\thefootnotetwo=\thefootnote%
	\renewcommand{\thefootnote}{}%
	\long\def\@makefntext##1{%
		\parindent 0em
		\noindent\hb@xt@0em{}##1%
	}%
	\ifthenelse{\equal{#1}{}}%
		{}%
		{\footnotetext{#1}}%
	\let\thefootnote=\thefootnotetwo%
	\mymakefntext%
}

\newcommand{\makethankshere}{\makethanks{\mbox{ }\\[-.75\baselineskip]\@researchsupported}}

%-------------------%

\newlength{\abovetitle}
\setlength{\abovetitle}{0.5\baselineskip}

\def\@maketitle{%
	\newpage%
	\thispagestyle{mytitlepage}%
	\vspace*{\abovetitle}%
	\let\footnote\thanks%
	{%
	  \centering%
      \ifthenelse{\boolean{jnetdy}}%		
      {\parbox[t]{\textwidth}{\centering\Large\bfseries\boldmath\@title}}%
      {\parbox[t]{.9\textwidth}{\centering\Large\bfseries\boldmath\@title}}%
      \par%
	  \vspace*{\baselineskip}
      \ifthenelse{\boolean{noauthor}}%
        %% For articles without any authors.
        %% Then \dedicated is intended to hold, e.g., venue and date of a conference.
        {\ifthenelse{\equal{\@dedicated}{}}{}%
          {\parbox[t]{.9\textwidth}{\centering\normalfont\normalsize\@dedicated\par}\par}%
        }%
        %% For research articles with authors.
        {\parbox[t]{.9\textwidth}{%
		 \centering\large%
         \ifthenelse{\boolean{prefaceetc}\or\equal{\@communicated}{}}%
           {\@author\par}%
           {\@author\\[2mm]{\small Communicated by \@communicated}\par}
         \ifthenelse{\equal{\@dedicated}{}}{}%
           {\vspace*{4mm}{\small\textit{\@dedicated}\par}}
       }\par%
       }
	}%
	\vspace*{2\baselineskip}
	%
	\ifthenelse{\boolean{jnetdy}\or\boolean{prefaceetc}}%
		{}%
		{%
			\noindent%
			{%
				\small\textbf{\abstractname.}~%
				\ifthenelse{\equal{\@abstract}{}}%
					{Supply an abstract using \string\abstract.}%
					{\@abstract}%
				\par%
			}%
			\vspace*{.7em}%
		}
	%
	\ifthenelse{\equal{\@keywords}{}}%
		{}%
		{\noindent{\small\raggedright\textbf{\keywordsname.}~\@keywords%
         \unskip.\par}\vspace*{.7em}}
	%
	\ifthenelse{\equal{\@classification}{}}%
		{}%
		{\noindent{\small\raggedright\textbf{\classificationname.}~\@classification%
         \unskip.\par}\vspace*{.7em}}
	%
	\vspace*{\baselineskip}
}% of \def\@maketitle

\newcommand\makeacknowledgments{
	\ifthenelse{\equal{\@acknowledgments}{}}%
		{}%
		{%
			\vspace{2.5ex \@plus .5ex \@minus .2ex}%
			\noindent\textbf{\acknowledgmentsname.}\hspace*{.5em}\xprint{\@acknowledgments}%
		}
}

\newcommand{\@receivedrevisedaccepted}{%
	Received~\xprint{\@received}%
	\ifthenelse{\equal{\@revised}{}}%
		{}%
		{;\hspace{2mm}revised~\xprint{\@revised}}% %% Don't delete the %-sign at the end of each line!
	%
	\ifthenelse{\equal{\@accepted}{}}%
		{}%
		{;\hspace{2mm}accepted~\xprint{\@accepted}}.%
	%
	\ifthenelse{\equal{\@publishedonline}{}}%
		{}%
		{\newline Published online~\xprint{\@publishedonline}.}%
}

\newcommand{\xprint}[1]{\ignorespaces #1\unskip}

\newcommand{\setaddress}[5][]{%
	\par\noindent%
	\begin{minipage}[b]{\textwidth}
		\ifthenelse{\equal{#1}{}}%
			{}%
			{\textbf{#1}\\[2mm]}%
		\raggedright%
	\xprint{#2}, \xprint{#3}, \xprint{#4}.
	%
	\ifthenelse{\equal{#5}{}}%
		{}%
		{\newline E-mail:~{\href{mailto:#5}{\tt #5}}}%
	%
	\end{minipage}\newline
	\addtocounter{numberauthors}{-1}
}

\newcommand{\makeauthordata}[1]{%
	\small%
    \vspace{3.5ex \@plus .5ex \@minus .2ex}%
    \ifthenelse{\boolean{nohistory}}%
      {}%
      {%
        \noindent\@receivedrevisedaccepted\par
        \vspace{3.5ex \@plus .5ex \@minus .2ex}%
      }%
	\setaddress[Author information]{\@authorone}{\@addressone}{\@countryone}{\@emailone}
	\ifthenelse{\equal{\value{numberauthors}}{0}}%
		{\label{#1}}%
		{%
			\setaddress{\@authortwo}{\@addresstwo}{\@countrytwo}{\@emailtwo}
			\ifthenelse{\equal{\value{numberauthors}}{0}}%
				{\label{#1}}%
				{%
					\setaddress{\@authorthree}{\@addressthree}{\@countrythree}{\@emailthree}
					\ifthenelse{\equal{\value{numberauthors}}{0}}%
					{\label{#1}}%
					{%
						\setaddress{\@authorfour}{\@addressfour}{\@countryfour}{\@emailfour}
						\ifthenelse{\equal{\value{numberauthors}}{0}}%
							{\label{#1}}%
							{\setaddress{\@authorfive}{\@addressfive}{\@countryfive}{\@emailfive}\label{#1}}
					}
				}
		}
    \enlargethispage{\baselineskip}  %NEU 02.07.2011
}

\def\makewarnings{
	\ifthenelse{\equal{\value{titletoolong}}{1}}%
		{%
			\typeout{}%
			\typeout{------------------------------------------------}%
			\typeout{WARNING: Title too long for running head.}%
			\typeout{PLEASE  supply a shorter form with \string\headlinetitle}
			\typeout{------------------------------------------------}%
			\typeout{}%
			\vspace{20mm} \par \large WARNING: Title too long for running head.
			\par \large PLEASE  supply a shorter form with \string\headlinetitle
		}%
		{}
	\ifthenelse{\equal{\value{authortoolong}}{1}}%
		{%
			\typeout{}%
			\typeout{------------------------------------------------}%
			\typeout{WARNING: Author names too long for running head.}%
			\typeout{PLEASE  supply a shorter form with \string\headlineauthor}
			\typeout{------------------------------------------------}%
			\typeout{}%
			\vspace{20mm} \par \large WARNING: Author names too long for running head.
			\par \large PLEASE  supply a shorter form with \string\headlineauthor
		}%
		{}
}% of \def\makewarnings

\let\origdocument\document
\let\endorigdocument\enddocument

\def\document{%
  \origdocument%
  \setcounter{page}{\value{startpage}}%
  \bibliographystyle{degruyter-plain}%
  %
  \ifthenelse{\boolean{noheadline}}%
   {}%
   {%
     \maketitle%
	 \ifthenelse{\boolean{nothanks}}{}{\makethanks{\@researchsupported}}%
   }
}

\ifthenelse{\boolean{noheadline}\or\boolean{noauthor}}%
  {\def\enddocument{\label{LastPage}\makewarnings\endorigdocument}}%
  {\def\enddocument{\makeauthordata{LastPage}\makewarnings\endorigdocument}}

\newcounter{section}
\newcounter{subsection}[section]
\newcounter{subsubsection}[subsection]
\newcounter{paragraph}[subsubsection]
\newcounter{subparagraph}[paragraph]

%-----------%
% Equations %
%-----------%

\renewcommand{\thesection}{\@arabic\c@section}
\renewcommand{\thesubsection}{\thesection.\@arabic\c@subsection}
\renewcommand{\thesubsubsection}{\thesubsection.\@arabic\c@subsubsection}

%% equation numbers like (1.1), (1.2), (1.3)
\newcommand{\SectionEquationNumbers}{
  \@addtoreset{equation}{section}%
  \renewcommand{\theequation}{\thesection.\arabic{equation}}%
}

%% equation numbers like (1), (2), (3)
\newcommand{\NoSectionEquationNumbers}{
  \renewcommand{\theequation}{\arabic{equation}}%
}

\ifthenelse{\boolean{seceqcount}}
	{\SectionEquationNumbers}
	{\NoSectionEquationNumbers}

\newcommand{\section}{%
	\@startsection{section}{1}{\z@}%
                 {-3.8ex \@plus -1ex \@minus -.2ex}%
                 {1.9ex \@plus.2ex}%
                 {\raggedright\normalfont\large\bfseries\boldmath}}
%
\newcommand{\subsection}{%
	\@startsection{subsection}{2}{\z@}%
                 {-3.5ex\@plus -1ex \@minus -.2ex}%
                 {1.4ex \@plus .2ex}%
                 {\raggedright\normalfont\normalsize\bfseries\boldmath}}
%
\newcommand{\subsubsection}{%
	\@startsection{subsubsection}{3}{\z@}%
                 {-3.3ex\@plus -1ex \@minus -.2ex}%
                 {1.3ex \@plus .2ex}%
                 {\raggedright\normalfont\normalsize\bfseries\boldmath}}
%
\newcommand{\paragraph}{%
	\@startsection{paragraph}{4}{\z@}%
                 {2.6ex \@plus1ex \@minus.2ex}% %rather than 3.25ex
                 {-1em}%
                 {\normalfont\normalsize\bfseries\boldmath}}
%
\newcommand{\subparagraph}{%
	\@startsection{subparagraph}{5}{\parindent}%
                 {3.25ex \@plus1ex \@minus .2ex}%
                 {-1em}%
                 {\normalfont\normalsize\bfseries\boldmath}}

\setlength{\leftmargini}{8.1mm} % de Gruyter magic number
\leftmargin\leftmargini
\setlength{\leftmarginii}{2.2em}
\setlength{\leftmarginiii}{1.87em}
\setlength{\leftmarginiv}{1.7em}
\setlength{\leftmarginv}{1em}
\setlength{\leftmarginvi}{1em}
\setlength{\labelsep}{.5em} % de Gruyter
\setlength{\labelwidth}{\leftmargini}
\addtolength{\labelwidth}{-\labelsep} % labelwidth + \labelsep = leftmargini

\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty

\renewcommand{\theenumi}{\@roman\c@enumi}
\renewcommand{\theenumii}{\@alph\c@enumii}
\renewcommand{\theenumiii}{\@roman\c@enumiii}
\renewcommand{\theenumiv}{\@Alph\c@enumiv}
\newcommand{\labelenumi}{\textup{(\theenumi)}}
\newcommand{\labelenumii}{\textup{\theenumii.}}
\newcommand{\labelenumiii}{\textup{\theenumiii.}}
\newcommand{\labelenumiv}{\textup{\theenumiv.}}
\renewcommand{\p@enumii}{\theenumi}
\renewcommand{\p@enumiii}{\theenumi(\theenumii)}
\renewcommand{\p@enumiv}{\p@enumiii\theenumiii}
%
\newcommand{\tinybullet}{{\tiny\raisebox{0.6mm}{$\bullet$}}}
%
\newcommand{\labelitemi}{\tinybullet}
\newcommand{\labelitemii}{\normalfont\bfseries\textendash}
\newcommand{\labelitemiii}{\textasteriskcentered}
\newcommand{\labelitemiv}{\textperiodcentered}
%
% Redefine the enumerate environment to get upright labels when using
% the ``optional argument'' construction provided by the enumerate package.
%
\let\origenumerate\enumerate
\let\endorigenumerate\endenumerate
\renewenvironment{enumerate}[1][no@ption]%
	{%
		\ifthenelse{\equal{#1}{no@ption}}%
		{\begin{origenumerate}}%
		{\begin{origenumerate}[\upshape#1]}%
	}%
	{\end{origenumerate}}
%
\newenvironment{xitemize}%
	{%   %% Unindented itemization.
		\begin{list}{\tinybullet}{%
			\setlength\labelwidth{0pt}
			\setlength\labelsep{7pt}
			\setlength\itemindent{7pt}
			\setlength\leftmargin{0pt}
			\setlength\listparindent{\parindent}
			\setlength\topsep{0.5\baselineskip}
			\setlength\partopsep{\parskip}
			\setlength\itemsep{0.5\baselineskip}
			\setlength\parsep{\parskip}
		}
	}%
	{\end{list}}
%
\newenvironment{yitemize}%   %% As xitemize, but with indention of the labels. Should be avoided.
	{%
		\begin{list}{\tinybullet}{%
			\setlength\labelwidth{-\parindent}
			\setlength\labelsep{7pt}
			\setlength\itemindent{7pt}
			\setlength\leftmargin{0pt}
			\setlength\listparindent{\parindent}
			\setlength\topsep{\parskip}
			\setlength\partopsep{\parskip}
			\setlength\itemsep{\parskip}
			\setlength\parsep{\parskip}
		}
	}%
	{\end{list}}
%
\newenvironment{xdescription}[1][(vii)]%   %% The parameter denotes the widest label.
	{%
		\begin{list}{\labelenumi}{
			\settowidth\labelwidth{\xdescriptionlabel{#1}}
			\setlength\leftmargin{\labelwidth+\labelsep}
			\let\makelabel\xdescriptionlabel
		}%
	}%
	{\end{list}}
%
\newcommand*{\xdescriptionlabel}[1]{\normalfont #1}
%
\newenvironment{ydescription}[1][(vii)]%   %% As xdescription, but with right-aligned labels.
	{%
		\begin{list}{\labelenumi}{
			\settowidth\labelwidth{\ydescriptionlabel{#1}}
			\setlength\leftmargin{\labelwidth+\labelsep}
			\let\makelabel\ydescriptionlabel
		}%
	}%
	{\end{list}}
%
\newcommand*{\ydescriptionlabel}[1]{\hfill\normalfont #1}

\newenvironment{description}%
	{%
		\list{}{%
			\labelwidth\z@ \itemindent-\leftmargin
			\let\makelabel\descriptionlabel
		}
	}%
	{\endlist}
%
\newcommand*{\descriptionlabel}[1]{\hspace\labelsep\normalfont\bfseries\boldmath #1}

% enumitem needs a predefind description enviroment.

%\usepackage{enumitem}

%\setlist{noitemsep}
%\setlist[1]{\labelindent=\parindent} % << Usually a good idea
%\setitemize{leftmargin=*}
%\setitemize[1]{label=$\triangleleft$}

%\setenumerate{labelsep=*, label=\textup{\Roman*}, leftmargin=1.5pc}
%\SetEnumerateShortLabel{i}{\textsc{\roman*}}
%\setenumerate[1]{label=\arabic*., ref=\arabic*}
%\setenumerate[2]{label=\emph{\alph*}), ref=\theenumi.\emph{\alph*}}
%\setenumerate[3]{label=\roman*), ref=\theenumii.\roman*}
%\setdescription{font=\sffamily\bfseries}

\newenvironment{verse}%
	{%
		\let\\\@centercr
		\begin{list}{}{
			\setlength{\itemsep}{\z@}%
			\setlength{\topsep}{1.5mm}%
			\setlength{\itemindent}{-1.5em}%
			\setlength{\listparindent}{\itemindent}%
			\setlength{\rightmargin}{\leftmargin}%
			\advance\leftmargin 1.5em%
		}%
		\item\relax%
	}
	{\end{list}}
%
\newenvironment{quotation}%
	{%
		\begin{list}{}{%
			\setlength{\listparindent}{1.5em}%
			\setlength{\topsep}{1.5mm}%
			\setlength{\itemindent}{\listparindent}%
			\setlength{\rightmargin}{\leftmargin}%
			\setlength{\parsep}{\z@\@plus\p@}%
		}%
		\item\relax}
	{\end{list}}
%
\newenvironment{quote}%
	{%
		\begin{list}{}{%
			\setlength{\topsep}{1.5mm}%
			\setlength{\rightmargin}{\leftmargin}%
		}%
		\item\relax%
	}
	{\end{list}}
%
\newenvironment{block}{\begin{trivlist}\item{}}{\end{trivlist}}

\newcommand{\appendix}{%
	\par%
	\setcounter{section}{0}%
	\setcounter{subsection}{0}%
	\gdef\thesection{\@Alph\c@section}%
}
%
\setlength{\arraycolsep}{5\p@}
\setlength{\tabcolsep}{6\p@}
\setlength{\arrayrulewidth}{0.4\p@}
\setlength{\doublerulesep}{2\p@}
\setlength{\tabbingsep}{\labelsep}
%
\skip\@mpfootins = \skip\footins
\setlength{\fboxsep}{3\p@}
\setlength{\fboxrule}{.4\p@}
\newcounter{figure}
\renewcommand{\thefigure}{\@arabic\c@figure}
\def\fps@figure{tb}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\upshape\figurename\nobreakspace\thefigure}
%
\newenvironment{figure}
	{\@float{figure}}
	{\end@float}
%
\newenvironment{figure*}
	{\@dblfloat{figure}}
	{\end@dblfloat}
%
\def\@floatboxreset{\global\@minipagefalse \centering}
%
\newcounter{table}
\renewcommand{\thetable}{\@arabic\c@table}
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename\nobreakspace\thetable}
%
\newenvironment{table}
	{\@float{table}}
	{\end@float}
%
\newenvironment{table*}
	{\@dblfloat{table}}
	{\end@dblfloat}
%
\let\l@table\l@figure

%----------------%
% Float captions %
%----------------%

\RequirePackage[style=base]{caption}
\captionsetup{margin=\parindent,justification=justified,%
  singlelinecheck=true,font=small,labelsep=period}
\usepackage[margin=0pt]{subfig}
% Although the caption package gives a compatibility warning,
% it is supported by our document class.

\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
%----------------------------------------------------------%
%          >> Simon Albroscheit am 16.01.2009 <<           %
%                                                          %
% Entfernt man in der Definition                           %
\DeclareRobustCommand*\cal{\@fontswitch\relax\mathcal}
% das ``\@fontswitch'', so liefert \cal sowohl mit als     %
% auch ohne mtpro2 kalligrafische Buchstaben!              %
%                                                          %
% Allerdings wird in $\cal A = A$ dann nur das erste A     %
% kalligrafisch, wohingegen mit der klassischen Definition %
% (und ohne mtpro2) BEIDE A's kalligrafisch werden --      %
% was natuerlich unschoen ist, aber verbreitet.            %
% Daher definieren wir:                                    %
%
\ifthenelse{\boolean{MTPro}}%
  {\renewcommand\cal{\@latex@error{%
   REPLACE \string\cal\ BY \mathcal}\mathcal}}%
  {}%
%
% Im Falle einer Warnung, sollte man also genau            %
% ueberpruefen, WELCHE Buchstaben der Autor kalligrafisch  %
% haben moechte und anschliessend in der Praeambel         %
% definieren: \renewcommand\cal{\mathcal}                  %
%----------------------------------------------------------%
\DeclareRobustCommand*\mit{\@fontswitch\relax\mathnormal}
%
\ifthenelse{\boolean{MTPro}}%
	{}%
	{%
		\DeclareSymbolFont{rmmath}{OT1}{\rmdefault}{m}{n}
		\SetSymbolFont{rmmath}{bold}{OT1}{\rmdefault}{b}{n}
		\DeclareMathSymbol{1}{\mathalpha}{rmmath}{"31}
		\DeclareMathSymbol{2}{\mathalpha}{rmmath}{"32}
		\DeclareMathSymbol{3}{\mathalpha}{rmmath}{"33}
		\DeclareMathSymbol{4}{\mathalpha}{rmmath}{"34}
		\DeclareMathSymbol{5}{\mathalpha}{rmmath}{"35}
		\DeclareMathSymbol{6}{\mathalpha}{rmmath}{"36}
		\DeclareMathSymbol{7}{\mathalpha}{rmmath}{"37}
		\DeclareMathSymbol{8}{\mathalpha}{rmmath}{"38}
		\DeclareMathSymbol{9}{\mathalpha}{rmmath}{"39}
		\DeclareMathSymbol{0}{\mathalpha}{rmmath}{"30}
		\SetSymbolFont{operators}{normal}{OT1}{ptmcm}{m}{n}
		%%\SetSymbolFont{operators}{bold}{OT1}{ptmcm}{bx}{n} %dann waeren grosse Griechen nicht fett in fetter Umgebung
	}

\DeclareSymbolFont{SFmathgreek}{OT1}{txss}{m}{n}   %TX Sans
\SetSymbolFont{SFmathgreek}{bold}{OT1}{txss}{b}{n} %TX Sans bold
\DeclareMathSymbol{\sfSigma}{\mathalpha}{SFmathgreek}{"06}
\DeclareMathSymbol{\sfGamma}{\mathalpha}{SFmathgreek}{"00}
\DeclareMathSymbol{\sfDelta}{\mathalpha}{SFmathgreek}{"01}
\DeclareMathSymbol{\sfTheta}{\mathalpha}{SFmathgreek}{"02}
\DeclareMathSymbol{\sfLambda}{\mathalpha}{SFmathgreek}{"03}
\DeclareMathSymbol{\sfXi}{\mathalpha}{SFmathgreek}{"04}
\DeclareMathSymbol{\sfPi}{\mathalpha}{SFmathgreek}{"05}
\DeclareMathSymbol{\sfUpsilon}{\mathalpha}{SFmathgreek}{"07}
\DeclareMathSymbol{\sfPhi}{\mathalpha}{SFmathgreek}{"08}
\DeclareMathSymbol{\sfPsi}{\mathalpha}{SFmathgreek}{"09}
\DeclareMathSymbol{\sfOmega}{\mathalpha}{SFmathgreek}{"0A}

\newcommand{\@pnumwidth}{1.55em}
\newcommand{\@tocrmarg}{2.55em}
\newcommand{\@dotsep}{4.5}
\setcounter{tocdepth}{2}

\newcommand{\tableofcontents}{\vspace{\baselineskip}{\small\@starttoc{toc}}}

\newcommand*{\l@section}[2]{%
	\ifnum \c@tocdepth >\m@ne
		\addpenalty{-\@highpenalty}%
		\vskip 0.4em \@plus\p@
		\setlength\@tempdima{1.5em}%
		\begingroup
		\parindent \z@ \rightskip \@pnumwidth
		\parfillskip -\@pnumwidth
		\leavevmode
		\advance\leftskip\@tempdima
		\hskip -\leftskip
		{\bfseries #1}\nobreak\
		\leaders\hbox{$\m@th\mkern \@dotsep mu\hbox{.}\mkern \@dotsep mu$}\hfil\nobreak\hb@xt@\@pnumwidth{\hss #2}\par
		\vskip 0.0em \@plus\p@
		\penalty\@highpenalty
		\endgroup
	\fi%
}

\newcommand*{\l@subsection}{\@dottedtocline{2}{1.5em}{2.3em}}
\newcommand*{\l@subsubsection}{\@dottedtocline{3}{3.8em}{3.2em}}
\newcommand*{\l@paragraph}{\@dottedtocline{4}{7.0em}{4.1em}}
\newcommand*{\l@subparagraph}{\@dottedtocline{5}{10em}{5em}}

%--------------------------------------------%
% for articles with class option ``noheadline'' %
%--------------------------------------------%
\newcommand{\mysection}[1]{%
	\vspace*{\baselineskip}
	{%
		\centering
		{\Large\textbf{#1}}\par
	}%
	\vspace*{4\baselineskip}%
}

\newcommand{\preface}{%
	\newpage
	\thispagestyle{empty}%
	\mysection{\prefacename}
	\markboth{\prefacename}{\prefacename}%
}

\newcommand{\xtableofcontents}{%
	\newpage
	\thispagestyle{empty}%
	\mysection{\contentsname}
	\markboth{\contentsname}{\contentsname}%
}

\newcommand{\mydotfill}{\cleaders\hbox{$\m@th\mkern3.9mu.\mkern3.9mu$}\hfill}

\newcommand\papertoc[3]{%
	\noindent\normalsize%
	\begin{minipage}[t]{\textwidth}
		\begin{minipage}[b]{.94\textwidth}
			\raggedright\textsc{#1}\\[.2ex]
			#2\phantom{I}\mydotfill\phantom{g}%
		\end{minipage}%
		\hfill#3%
	\end{minipage}%
\vspace{3mm}%
}

\newcommand{\othertoc}[2]{%
	\noindent\normalsize%
	\begin{minipage}[t]{\textwidth}
		\begin{minipage}[b]{.94\textwidth}
			\raggedright%
			#1\phantom{I}\mydotfill\phantom{g}
		\end{minipage}%
	\hfill#2%
	\end{minipage}%
	\vspace{3mm}%
}

\newenvironment{theindex}%
	{
		\newpage
		\thispagestyle{empty}
		\if@twocolumn
		\@restonecolfalse
		\else
		\@restonecoltrue
		\fi
		\columnseprule \z@
		\columnsep 35\p@
		\twocolumn[\vspace*{\baselineskip}\mysection{\indexname}]
		\markboth{\indexname}{\indexname}
		\parindent\z@
		\parskip\z@ \@plus .3\p@\relax
		\let\item\@idxitem\small
	}%
	{%
		\if@restonecol\onecolumn
		\else\clearpage\fi%
	}

%--------------------------------------------%

\newdimen{\bibindent}
\setlength{\bibindent}{1.5em}

\newenvironment{thebibliography}[1]%
	{
		\makeacknowledgments
		\section*{\refname}%
		\addcontentsline{toc}{section}{\refname}
		\small
		\list{\@biblabel{\@arabic\c@enumiv}}{
			\settowidth\labelwidth{\@biblabel{#1}}%
			\leftmargin\labelwidth
			\advance\leftmargin\labelsep
			\@openbib@code
			\usecounter{enumiv}%
			\let\p@enumiv\@empty
			\renewcommand\theenumiv{\@arabic\c@enumiv}%
		}%
		%
		\sloppy
		\clubpenalty4000
		\@clubpenalty \clubpenalty
		\widowpenalty4000%
		\sfcode`\.\@m%
	}
	{%
		\def\@noitemerr{\@latex@warning{Empty `thebibliography' environment}}%
		\endlist%
	}

\newcommand{\newblock}{\hskip .11em\@plus.33em\@minus.07em}
\let\@openbib@code\@empty

\newcommand{\@idxitem}{\par\hangindent 40\p@}
\newcommand{\subitem}{\@idxitem \hspace*{20\p@}}
\newcommand{\subsubitem}{\@idxitem \hspace*{30\p@}}
\newcommand{\indexspace}{\par \vskip 10\p@ \@plus5\p@ \@minus3\p@\relax}

\renewcommand\footnoterule{%
	\kern-3\p@
	\hrule\@width.4\columnwidth
	\kern2.6\p@%
}

\newcommand{\abstractname}{Abstract}
\newcommand{\acknowledgmentsname}{Acknowledgments}
\newcommand{\appendixname}{Appendix}
\newcommand{\bibname}{\refname}

\newcommand{\contentsname}{Contents}
\newcommand{\figurename}{Figure}
\newcommand{\indexname}{Index}
\newcommand{\keywordsname}{Keywords}
\newcommand{\partname}{Part}
\newcommand{\prefacename}{Preface}
\newcommand{\refname}{Bibliography}
\newcommand{\sectionname}{Section}
\newcommand{\tablename}{Table}
%
\def\today{\ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space\number\day, \number\year}

\pagenumbering{arabic}

\if@twoside
\else
	\raggedbottom
\fi

\if@twocolumn
	\twocolumn
	\sloppy
	\flushbottom
\else
	\onecolumn
\fi


%---------------------------------%
% Changes for the amsthm-package. %
%---------------------------------%

\renewcommand{\thmhead@plain}[3]{%
	\thmname{#1}%
	\thmnumber{\@ifnotempty{#1}{ }\@upn{#2}}%
	\thmnote{ {\the\thm@notefont(#3)\the\thm@headpunct}}% NEW
}

\let\thmhead\thmhead@plain

\def\swappedhead#1#2#3{%
	\thmnumber{#2}%
	\thmname{\@ifnotempty{#2}{~}#1}%
	\thmnote{ {\the\thm@notefont(#3)\the\thm@headpunct}}% NEW
}

\let\swappedhead@plain=\swappedhead

\def\@begintheorem#1#2[#3]{%
	\deferred@thm@head{%
		\the\thm@headfont \thm@indent
		\@ifempty{#1}{\let\thmname\@gobble}{\let\thmname\@iden}%
		\@ifempty{#2}{\let\thmnumber\@gobble}{\let\thmnumber\@iden}%
		\@ifempty{#3}{\let\thmnote\@gobble}{\let\thmnote\@iden}%
		\thm@swap\swappedhead\thmhead{#1}{#2}{#3}%
		\@ifempty{#3}{\the\thm@headpunct}{}%        %NEW
		%\the\thm@headpunct                         %NEW
		\thmheadnl%
		\hskip\thm@headsep
	}%
	\ignorespaces%
}

\newlength\qedrulewidth
\setlength\qedrulewidth{.37pt} %% magic
\newlength\qedwidth
\setlength\qedwidth{1ex} %% magic
\newlength\qedheight
\setlength\qedheight{\qedwidth-2\qedrulewidth}
\newlength\qedinnerwidth
\setlength\qedinnerwidth{\qedwidth-2\qedrulewidth}
%
\renewcommand{\qedsymbol}{%
	\vbox{%
		\hrule width\qedwidth height\qedrulewidth% %upper bar
		\hbox{%
			\vrule height\qedheight width\qedrulewidth% %left bar
			\hskip\qedinnerwidth%
			\vrule width\qedrulewidth% right bar
		}%
		\hrule width\qedwidth height\qedrulewidth% lower bar
	}
}

\newlength{\preskip}
\setlength{\preskip}{12\p@ \@plus5\p@ \@minus3\p@}

\newlength{\postskip}
\setlength{\postskip}{8\p@ \@plus3\p@ \@minus1.5\p@}

\newlength{\postskipp}
\setlength{\postskipp}{11\p@ \@plus3\p@ \@minus1.5\p@}

\newtheoremstyle{plain}%
	{\preskip}{\postskip}  %Abstand oben und unten
	{\itshape}             %Textlayout
	{0pt}                  %Einzug
	{\bfseries}{.}{ }      %"Uberschrift und Abstand
	{}                     %=Formatierungen der "Uberschrift wie bei plain

\newtheoremstyle{definition}%
	{\preskip}{\postskipp}  %Abstand oben und unten
	{\normalfont}          %Textlayout
	{0pt}                  %Einzug
	{\bfseries}{.}{ }      %"Uberschrift und Abstand
	{}                     %=Formatierungen der "Uberschrift wie bei plain


\ifthenelse{\boolean{integers}}%
	{\theoremstyle{plain}\newtheorem{theorem}{Theorem}}%
	{\theoremstyle{plain}\newtheorem{theorem}{Theorem}[section]}

\theoremstyle{plain}
	\newtheorem*{theorem*}{Theorem}
	\newtheorem{claim}[theorem]{Claim}
	\newtheorem*{claim*}{Claim}
	\newtheorem{corollary}[theorem]{Corollary}
	\newtheorem*{corollary*}{Corollary}
	\newtheorem{lemma}[theorem]{Lemma}
	\newtheorem*{lemma*}{Lemma}
	\newtheorem{proposition}[theorem]{Proposition}
	\newtheorem*{proposition*}{Proposition}
\theoremstyle{definition}
	%\newtheorem{algorithm}[theorem]{Algorithm} % not compatible with algorithm2e
	%\newtheorem*{algorithm*}{Algorithm}
	\newtheorem{assumption}[theorem]{Assumption}
	\newtheorem*{assumption*}{Assumption}
	\newtheorem{conjecture}[theorem]{Conjecture}
	\newtheorem*{conjecture*}{Conjecture}
	\newtheorem{definition}[theorem]{Definition}
	\newtheorem*{definition*}{Definition}
	\newtheorem{example}[theorem]{Example}
	\newtheorem*{example*}{Example}
	\newtheorem{hypothesis}[theorem]{Hypothesis}
	\newtheorem*{hypothesis*}{Hypothesis}
    \newtheorem{notation}[theorem]{Notation}
    \newtheorem*{notation*}{Notation}
	\newtheorem{problem}[theorem]{Problem}
	\newtheorem*{problem*}{Problem}
	\newtheorem{property}[theorem]{Property}
	\newtheorem*{property*}{Property}
	\newtheorem{remark}[theorem]{Remark}
	\newtheorem*{remark*}{Remark}
%---------------------------------%

\newcounter{mynote}
\setcounter{mynote}{0}
\newcommand{\mynote}[2][\themynote]{%
	\refstepcounter{mynote}%
	\label{mynote:#1}%
	\mbox{}\marginpar{%
		\hspace{0pt}%
		\RaggedRight{%
			\sffamily\tiny\bfseries Note~\themynote:%
		}\\ %
	\tiny#2}%
}

\endinput
%%
%% End of the file `degruyter-journal-a.cls'.
