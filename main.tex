% Copyright 2004 by Till Tantau <tantau@users.sourceforge.net>.
%
% In principle, this file can be redistributed and/or modified under
% the terms of the GNU Public License, version 2.
%
% However, this file is supposed to be a template to be modified
% for your own needs. For this reason, if you use this file as a
% template and not specifically distribute it as part of a another
% package/program, I grant the extra permission to freely copy and
% modify this file as you see fit and even to delete this copyright
% notice. 

\documentclass[aspectratio=169, handout]{beamer}

% There are many different themes available for Beamer. A comprehensive
% list with examples is given here:
% http://deic.uab.es/~iblanes/beamer_gallery/index_by_theme.html
% You can uncomment the themes below if you would like to use a different
% one:
%\usetheme{AnnArbor}
%\usetheme{Antibes}
%\usetheme{Bergen}
%\usetheme{Berkeley}
%\usetheme{Berlin}
\usetheme{Boadilla}
%\usetheme{boxes}
%\usetheme{CambridgeUS}
%\usetheme{Copenhagen}
%\usetheme{Darmstadt}
%\usetheme{default}
%\usetheme{Frankfurt}
%\usetheme{Goettingen}
%\usetheme{Hannover}
%\usetheme{Ilmenau}
%\usetheme{JuanLesPins}
%\usetheme{Luebeck}
% \usetheme{Madrid}
%\usetheme{Malmoe}
%\usetheme{Marburg}
%\usetheme{Montpellier}
%\usetheme{PaloAlto}
%\usetheme{Pittsburgh}
%\usetheme{Rochester}
%\usetheme{Singapore}
%\usetheme{Szeged}
%\usetheme{Warsaw}

\title[Coded Computing]{Coded Computing for Straggler Mitigation, Security and Privacy}
% A subtitle is optional and this may be deleted
\subtitle{EE605 Error Correcting Codes}

% \author{F.~Author\inst{1} \and S.~Another\inst{2}}
\author[Param, Anupam]{Param Rathour 190070049\\Anupam Nayak 19D070010}
% - Give the names in the same order as the appear in the paper.
% - Use the \inst{?} command only if the authors have different
%   affiliation.

\institute[IIT Bombay]{Department of Electrical Engineering\\
Indian Institue of Technology Bombay} % (optional, but mostly needed)
% {
%   \inst{1}%
%   Department of Computer Science\\
%   University of Somewhere
%   \and
%   \inst{2}%
%   Department of Theoretical Philosophy\\
%   University of Elsewhere}
% - Use the \inst command only if there are several affiliations.
% - Keep it simple, no one is interested in your street address.

\date{Autumn 2021-22}
% - Either use conference name or its abbreviation.
% - Not really informative to the audience, more for people (including
%   yourself) who are reading the slides online

\subject{Coding Theory}
% This is only inserted into the PDF information catalog. Can be left
% out. 

% If you have a file called "university-logo-filename.xxx", where xxx
% is a graphic format that can be processed by latex or pdflatex,
% resp., then you can add a logo as follows:

% \pgfdeclareimage[height=0.5cm]{university-logo}{university-logo-filename}
% \logo{\pgfuseimage{university-logo}}

% Delete this, if you do not want the table of contents to pop up at
% the beginning of each subsection:
% \AtBeginSubsection[]
% {
%   \begin{frame}<beamer>{Outline}
%     \tableofcontents[currentsection,currentsubsection]
%   \end{frame}
% }
\newcommand{\bF}{\mathbb{F}}
\newcommand{\bV}{\mathbb{V}}
\newcommand{\bU}{\mathbb{U}}
\newcommand{\bW}{\mathbb{W}}
\newcommand{\bR}{\mathbb{R}}
\newcommand{\floor}[1]{\lfloor #1 \rfloor}
\usepackage{ragged2e}
% \usepackage{etoolbox}
% \apptocmd{\frame}{}{\justifying}{} % Allow optional arguments after frame.
\renewcommand{\raggedright}{\leftskip=0pt \rightskip=0pt plus 0cm}
\apptocmd{\frame}{}{\justifying}{}
% \addtobeamertemplate{}{}{\justifying}
\setbeamersize{text margin left=2em,text margin right=3em}
% \beamerdefaultoverlayspecification{<+->}
% \addtobeamertemplate{proof begin}{%
%     \setbeamercolor{block title}{fg=red!50!black,bg=red!25!white}%
%     \setbeamercolor{block body}{fg=black, bg=red!10!white}%
% }{}
% Let's get started
\begin{document}

\begin{frame}
  \titlepage
\end{frame}

\begin{frame}{Outline}
  \tableofcontents
  % You might wish to add the option [pausesections]
\end{frame}

% Section and subsections will appear in the presentation overview
% and table of contents.
\section{Straggler Mitigation in Distributed Matrix Multiplications}

\subsection{Problem Formulation}

% \begin{frame}{First Slide Title}{Optional Subtitle}
%   \begin{itemize}
%   \pause\item {
%     My first point.
%   }
%   \pause\item {
%     My second point.
%   }
%   \end{itemize}
% \end{frame}


% You can reveal the parts of a slide one at a time
% with the \pause command:
% \begin{frame}{Second Slide Title}
%   \begin{itemize}
%   \pause\item {
%     First item.
%     \pause % The slide will pause after showing the first item
%   }
%   \pause\item {   
%     Second item.
%   }
%   % You can also specify when the content should appear
%   % by using <n->:
%   \pause\item<3-> {
%     Third item.
%   }
%   \pause\item<4-> {
%     Fourth item.
%   }
%   % or you can use the \uncover command to reveal general
%   % content (not just \pause\items):
%   \pause\item<5-> {
%     Fifth item. \uncover<6->{Extra text in the fifth item.}
%   }
%   \end{itemize}
% \end{frame}
\begin{frame}{Distributed Matrix Multiplication}{Problem formulation}
We have two input matrices
    \begin{itemize}
        \pause\item $A \in \mathbb{F}_q^{s \times r}$
        \pause\item $B \in \mathbb{F}_q^{s \times t}$ for some sufficiently large finite field $\mathbb{F}_q$
        \pause\item \textbf{To compute} $C = A^TB$. It is implicitly assumed in general that one of these matrices is tall.
    \end{itemize}
\pause Each worker has to be assigned fraction of the coded a fraction of the submatrix
    \begin{itemize}
    \pause \item Each of the N workers stores $\frac{1}{m}$ fraction of $A$ and $\frac{1}{n}$ fraction of $B$ where $m,n \in \mathbb{N}$
    \pause \item Thus $A_i \in \mathbb{F}_q^{s \times \frac{r}{m}}$ and $B_i \in \mathbb{F}_q^{s\times \frac{t}{n}}$
    \end{itemize}
Idea is for the master to use something like an MDS encoding to generate each of these sub matrices such that it has to wait only for $k$ of these workers(fastest) to generate the output in order to uniquely identify the product.
\end{frame}

\begin{frame}{Problem Formulation}{Key Ideas}
\begin{figure}
    \centering
    \includegraphics[width = 0.9\linewidth]{sys_tall_v2.pdf}
    \caption{Overview of the distributed matrix multiplication problem\footnote{\tiny{Q. Yu, M. A. Maddah-Ali, and A. S. Avestimehr, “Straggler mitigation in distributed matrix multiplication: Fundamental limits and optimal coding,” 2020.}}.}
    \label{fig:E2}
\end{figure}
\end{frame}
\begin{frame}{Computation Strategy}

    \begin{block}{Choice of functions}
    A computation strategy is defined as a set of $2N$ functions
    \[f = (f_0,f_1.....f_{N-1})\]
    \[g = (g_0,g_1.....g_{N-1})\]
   that are used to compute $A_i = f_i(A)$, $B_i = g_i(B)$    $\forall i \in \{0,1,2,...N-1\}$
    \end{block}
    \pause For any integer k we say that the system is $k$ recoverable if the
master can recover the product $C$ using output from any $k$ workers\pause\\
We define $k(f,g)$ as the least integer $k$ for which the system defined by $f,g$ is $k$ recoverable


\end{frame}
\subsection{Computation Strategy}
\begin{frame}{Computation Strategy}

 \begin{block}{Optimum Recovery Threshold}
    The lowest among the recovery thresholds across
all computation strategies
 \[K^* = \min_{f,g} k(f,g)\]
    \end{block}
\pause
State-of-the-art schemes include 
\begin{itemize}
    \item 1D MDS scheme
     \[K_{1D \;MDS} = N-\frac{N}{n}+m\]
     \item In the same paper an alternative scheme for the special case of $m =n$ referred to as the product code achieves a threshold of
     \[K_{prod} = 2(m-1)\sqrt{N} - (m-1)^2 +1 \]
\end{itemize}
\end{frame}

\begin{frame}{Main Result}
\begin{theorem}
The distributed matrix multiplication problem of computing $A^TB$ described above has a minimum recovery threshold of 
\[K^* = mn\]
\pause
Further $\exists$ a computation strategy referred to as polynomial code which achieves the above $K^*$ which allows for efficient decoding at the master node with computational complexity of polynomial decoding with $mn$ points.
    \end{theorem}
\pause
Decoding polynomials codes essentially a polynomial interpolation problem, which can be solved in time almost linear to the input size.This is enabled by designing the computing strategies such that the computed products form a Reed-Solomon code.
\end{frame}


% \begin{frame}{Comparison of the three thresholds}
% \begin{figure}
%     \centering
%     \includegraphics[width = 0.6\linewidth]{E1.png}
%     \caption{}
%     \label{fig:E1}
% \end{figure}
% \end{frame}

\begin{frame}{Comparison of the four thresholds}
\begin{figure}
    \centering
    \includegraphics[width = 0.5\linewidth]{compare_v6.pdf}
    \caption{Comparison\footnote{\tiny{Q. Yu, M. A. Maddah-Ali, and A. S. Avestimehr, “Straggler mitigation in distributed matrix multiplication: Fundamental limits and optimal coding,” 2020.}}}
    \label{fig:com}
\end{figure}
\end{frame}

\begin{frame}{Proof I}
A sketch of the proof
\vspace{3pt}\\
Given parameters $\alpha,\beta \in \mathbb{N} $, we define the $(\alpha,\beta)$-
polynomial code $\forall i \in \{0,1....N-1\}$ as
\uncover<2->{\[\overline{A_i} = \sum_{j=0}^{m-1}A_j x_i^{j\alpha}\]
\[\overline{B_i} = \sum_{j=0}^{n-1}B_j x_i^{j\beta}\]
\[\overline{C_i} = \overline{A_i}^T\overline{B_i} = \sum_{j=0}^{m-1}\sum_{k=0}^{n-1}A_j^TB_kx_i^{j\alpha+k\beta}\]}
\uncover<3->{Note here we need to carefully choose $\alpha$ and $\beta$ such that no two terms have the same power of $x$. One such choice being $\alpha = 1$, $\beta = m$. we define $h(x)$ as follows}
\end{frame}

\begin{frame}{Proof II}
\[h(x) =  \sum_{j=0}^{m-1}\sum_{k=0}^{n-1}A_j^TB_kx^{j+k\beta}\]
\pause
This is polynomial of degree $mn-1$. Since $x_i$ are chosen to be different in order to recover $C$ we need the output from any $mn$ workers which is essentially interpolating using $mn$ points, For even less complex decoding we may use the Reed Solomon decoding algorithm.
\vspace{3pt}\\
\pause So far we have a scheme that achieves the bound in the theorem. We need to prove that we require output from atleast $mn$ workers to recover $C$ in order to prove optimality



\end{frame}

\begin{frame}{Proof III}


A sketch of this proof is as follows
\vspace{3pt}\\
\begin{itemize}
  

  \item WLOG Let $A$ be an arbitrary fixed tall matrix
$(s\geq r)$ and $B$ is sampled from a uniform distribution over $\mathbb{F}_q^{s\times t}$
\uncover<2->{
  \item Thus one can check the distribution of $C = A^TB$ will be uniform over $\mathbb{F}_q^{r\times t}$}
\uncover<3->{  \item This means we need to recover a random variable with entropy 
$H(C) = rt \log q$
}
\uncover<4->{  \item Since each worker outputs $\frac{rt}{mn}$ elements of $\mathbb{F}_q$ it provides atmost $\frac{rt}{mn} \log q$ bits of information
\vspace{3pt}\\
hence $K>mn$
}

\end{itemize}
\end{frame}


\begin{frame}{Performance of the polynomial code on other evaluation metrics}
\begin{itemize}
    \item \textbf{Computation
latency} is defined as the amount of time required for the master to
collect enough information to decode C, For any other computation strategy we have
\[T \geq T_{poly}\]
\uncover<2->{\item \textbf{Probability of failure given a deadline} is defined as the probability
that the master does not receive enough information to decode C at
a predefined time t
\[P(T > t) \geq P(T_{poly} > t)\]}
\uncover<3->{\item \textbf{Communication load} is defined as the minimum number of bits needed to be extracted in
order to complete the computation. The below mentioned bound is achieved by the polynomial code.
\[L^* = rt \log_2 q\]}
\end{itemize}


\end{frame}


\section{Lagrange Coded Computing -- Going beyond Matrix Algebra}
\begin{frame}{Polynomial Evaluation}{Problem Formulation}
We have a distributed computing environment with a master and~$N$ workers
    \begin{itemize}
        \pause\item \textbf{Dataset} $X = (X_1,\ldots,X_k)$ where $X_i$ is a element of a vector space $\bV$ over $\bF$
        \pause\item $f: \bV \rightarrow \bU $ is a multivariate polynomial with vector coefficients and degree $=\deg f$
        \pause\item \textbf{To compute} $Y_1\triangleq f(X_1), \ldots, Y_K\triangleq f(X_K)$
    \end{itemize}
Each worker has already stored a fraction of the coded dataset prior to computation
    \begin{itemize}
    \pause\item The $i$\textsuperscript{th} worker stores $\tilde{X}_i\triangleq g_i(X_1,\ldots,X_K)$, where~$g_i$ is a (possibly random) function, refered to as the encoding function of that worker. (~$i\in [N]$ and ~$[N]\triangleq\{1,\ldots,N\}$)
    \pause\item Each worker~$i \in[N]$ computes~$\tilde{Y}_i\triangleq f(\tilde{X}_i)$ and returns the result to the master.
    \end{itemize}
The master waits for {a subset of fastest workers} and then decodes $Y_1,\ldots,Y_K$.

Linear encoding strategy gives simple yet concrete implmentation,
\end{frame}
\subsection{Polynomial Evaluations}
\begin{frame}{Some Polynomial Evaluations tasks}
    \begin{example}[Linear Computation]
    \textbullet\ Compute $A\vec{b}$ for some dataset $A=\{A_i\}_{i=1}^K$ and vector~$\vec{b}$
    \pause
    
    Let ~$\bV$ be the space of matrices over~$\bF$, $\bU$ be the space of vectors over~$\bF$, $X_i$ be $A_i$, and~$f(X_i)=X_i\cdot \vec{b}$ for all~$i\in[K]$. (suitable dimensions to be assigned to $\bV, \bU$)
    \end{example}
    \pause
    \begin{example}[Bilinear Computation]
    \textbullet\ Compute element-wise products $\{A_i\cdot B_i\}_{i=1}^{K}$ of two matrices $\{A_i\}_{i=1}^K$ and~$\{B_i\}_{i=1}^K$.
    \pause
    Let ~$\bV$ be the space of pairs of two matrices, $\bU$ be the space of matrices, $X_i=(A_i, B_i)$, and~$f(X_i)=A_i\cdot B_i$ for all~$i\in[K]$. (suitable dimensions to be assigned to $\bV, \bU$) 
    \end{example}
\end{frame}
\begin{frame}{Result on Minimum Workers required for recovery}
    \pause\begin{theorem}\label{thm:lccs}
        Given a number of workers~$N$ and a dataset~$X=(X_1,\ldots,X_K)$,  for distributedly computing $f$, the minimum recovery threshold is given by
        \begin{align}\label{eq:lccsThm}
            K^* = \begin{cases} 
        (K-1)\deg f+1 & K \deg f - 1 \leq N\\
        N - \floor{N/K} + 1& else
        \end{cases}
        \end{align} 
    \end{theorem}
    \pause\begin{theorem}
        The Lagrange Coded Computing is optimal i.e. it minimizes the recovery threshold
    \end{theorem}
    \begin{proof}
    Coming up.
    
    \end{proof}
\end{frame}
\begin{frame}{Polynomial Interpolation}\pause
    Given a set of $N + 1$ data points $(x_i, y_i)$ where no two $x_i$ are the same, a polynomial ${\displaystyle p:\mathbb {R} \rightarrow \mathbb {R} }$ is said to interpolate the data if ${\displaystyle p(x_{j})=y_{j}}$ for each ${\displaystyle j\in [N+1]}$\pause
    \begin{block}{Construction using System of Linear Equations}
    Suppose that the interpolation polynomial is in the form
    \begin{equation*}
        p(x)=a_{n}x^{n}+a_{n-1}x^{n-1}+\cdots +a_{2}x^{2}+a_{1}x+a_{0}
    \end{equation*}\pause
    Now, we can frame a system of linear equations as $p(x_i) = y_i \quad\mbox{for } i \in [N+1]$
    \begin{equation*}
        \begin{bmatrix}
x_0^n  & x_0^{n-1} & x_0^{n-2} & \ldots & x_0 & 1 \\
x_1^n  & x_1^{n-1} & x_1^{n-2} & \ldots & x_1 & 1 \\
\vdots & \vdots    & \vdots    &        & \vdots & \vdots \\
x_n^n  & x_n^{n-1} & x_n^{n-2} & \ldots & x_n & 1
\end{bmatrix}
\begin{bmatrix} a_n \\ a_{n-1} \\ \vdots \\ a_0 \end{bmatrix}  =
\begin{bmatrix} y_0 \\ y_1 \\ \vdots \\ y_n \end{bmatrix}.
    \end{equation*}
    \end{block}
\end{frame}
\begin{frame}{Lagrange Interpolation\footnote{\tiny{\url{https://commons.wikimedia.org/wiki/File:Lagrange_polynomial.svg}}}}
\begin{equation*}
    p(x)={\frac {(x-x_{1})(x-x_{2})\cdots (x-x_{n})}{(x_{0}-x_{1})(x_{0}-x_{2})\cdots (x_{0}-x_{n})}}y_{0}+\cdots +{\frac {(x-x_{0})(x-x_{1})\cdots (x-x_{n-1})}{(x_{n}-x_{0})(x_{n}-x_{1})\cdots (x_{n}-x_{n-1})}}y_{n}
\end{equation*}\pause
% \begin{figure}
%     \centering
%     \includegraphics[width = 0.5\linewidth]{Lagrange_polynomial.pdf}
%     % \caption{}
%     \label{fig:lp}
% \end{figure}
\begin{equation}
p(x) = \sum _{i=0}^{n}{\Bigg (}\prod _{\stackrel {\!0\,\leq \,j\,\leq \,n}{j\,\neq \,i}}{\frac {x-x_{j}}{x_{i}-x_{j}}}{\Bigg )}y_{i}
\qquad\qquad
    \vcenter{\hbox{\includegraphics[width=0.45\linewidth]{Lagrange_polynomial.pdf}}}
\end{equation}
\end{frame}
\begin{frame}{Lagrange Coded Computing}{Key Ideas}
    \begin{itemize}
        \pause\item The Lagrange interpolation polynomial is used to create encoding of the input dataset inserting computational redundancy in a coded form across the workers.
        \pause\item The computations at the each worker amount to evaluations of a composition of this polynomial with the desired function $f$ resulting in another polynomial $h_i$.
        \pause\item Decode $Y_1,\ldots, Y_K$ using only $K^*$ of $h_i$'s by evaluating each at certain points.
    \end{itemize}
\end{frame}
\begin{frame}{Lagrange Coded Computing}{Encoding}
    \begin{itemize}
        \pause\item Select any $K$ distinct elements $\beta_1,\ldots,\beta_K$ from $\bF$, and find a polynomial $u: \bF \rightarrow \bV$ of degree $K - 1$ such that $u(\beta_i) = X_i$ for $i \in [K]$.
        
        This can be accomplished by Lagrange interpolation polynomial
        \begin{equation}
            u(z) \triangleq \sum_{j\in[K]}X_j\cdot \prod_{k\in [K]\setminus\{j\}}\frac{z-\beta_k}{\beta_j-\beta_k}
        \end{equation}
        \pause\item Now, select $N$ distinct elements $\alpha_1,\ldots,\alpha_N$ from $\bF$ and  encode the input variables by letting $\tilde{X} = u(\alpha_i)$ for $i \in [N]$
        \begin{equation}
            \tilde{X}_i = g_i(X) = u(\alpha_i) \triangleq \sum_{j\in[K]}X_j\cdot \prod_{k\in [K]\setminus\{j\}}\frac{\alpha_i-\beta_k}{\beta_j-\beta_k}
        \end{equation}

    \end{itemize}
\end{frame}
\begin{frame}{Lagrange Coded Computing}{Decoding}
    \begin{itemize}
        \pause\item Each worker $i$ computes $\tilde{Y}_i = f(\tilde{X}_i) = f(u(\alpha_i))$ and sends $\tilde{Y}_i$ to the master.
        \pause\item This composition $f(u(z))$ is also a polynomial with degree $\leq (K-1)\deg f$.
        \pause\item Now, any $(K-1)\deg f+1$ workers return the evaluations at $(K-1)\deg f+1$ points. This gives a unique $f(u(z))$ which can be interpolated using Lagrange polynomials.
        \pause\item Then, the master evaluates it at~$\beta_i$ for every~$i\in[K]$ to obtain~$f(u(\beta_i))=f(X_i)$, %\blue
        % \pause\item Then, the master calculates $Y_i = f(X_i) = f(u(\beta_i))$ for $i \in [K]$ by evaluating $f(u(z))$.
    \end{itemize}
    \pause
    Note that if, number of workers are small ($N < K \deg f -1$), $K^*$ can be easily achieved by replicating every $X_i$ by atleast $\floor{N/K}$ times. Now, every set of $N-\floor{N/K}+1$ computation contains at least one copy of $f(X_i)$ for every $i$.
\end{frame}
\subsection{Secure and Private Multiparty Computing}
\begin{frame}{Secure and Private Multiparty Computing}{Overview}
    \begin{figure}
        \centering
        \includegraphics[width = 0.5\linewidth]{sys_v7.pdf}
        \caption{Overview of Secure and Private Multiparty Computing\footnote{\tiny J. S. Qian Yu, Netanel Raviv and A. S. Avestimehr, “Lagrange coded computing: Optimal design for resiliency, security and privacy,” 2018 }}
        \label{fig:sat}
    \end{figure}
\end{frame}
\begin{frame}{Secure and Private Multiparty Computing}{Terminology}
    \begin{itemize}
        \pause\item  \textbf{Resiliency} (robustness against stragglers) In a \textit{$S$-resilient} system, the master must be able to obtain the correct values of~$Y_1,\ldots,Y_K$ even if up to~$S$ workers delay/fail.
        \pause\item  \textbf{Security} (robustness against adversaries) In a \textit{$A$-secure} system, the master must be able to obtain correct values of~$Y_1,\ldots,Y_K$ even if up to~$A$ workers return arbitrarily erroneous results.
        \pause\item \textbf{Privacy} (robustness against collusion) In a \textit{$T$-private} system, the workers cannot infer anything about the content of the dataset, even if up to~$T$ of them collude, 
        
        Formally, for every $\mathcal{T} \subseteq [N]$ of size at most $T$, we must have $I(X;\tilde{X}_\mathcal{T})=0$
    \end{itemize}
    \pause The tuple $(S, A, T)$ is achievable if there exists an encoding and decoding scheme that can complete the computations in the presence of up to $S$ stragglers, up to $A$ adversarial workers, whilst keeping the dataset private against sets of up to $T$ colluding workers.
\end{frame}
\begin{frame}{Main Result}
The below theorem characterizes the region for $(S, A, T)$ that LCC achieves
    \pause\begin{theorem}\label{thm:lcc}
        Given a number of workers~$N$ and a dataset~$X=(X_1,\ldots,X_K)$, LCC provides an~$S$-resilient, $A$-secure, and~$T$-private scheme for computing~$\{f(X_i)\}_{i=1}^K$ for any polynomial $f$, as long as
        \begin{align}\label{eq:lccThm}
            (K+T-1)\deg f+S+2A+1\leq N.
        \end{align}  
    \end{theorem}
    \pause\begin{block}{Interesting Fact}
    One additional worker can increase its resiliency to stragglers by 1, or increase its robustness to adversaries by 1/2, while maintaining the privacy constraint. Sounds familiar?
    \end{block} 
\end{frame}
\begin{frame}{Lagrange Coded Computing}{Encoding}
    \begin{itemize}
        \pause\item Select any $K+T$ distinct elements $\beta_1,\ldots,\beta_{K+T}$ from $\bF$, and find a polynomial $u: \bF \rightarrow \bV$ of degree $K+T - 1$ such that $u(\beta_i) = X_i$ for $i \in [K]$ and $u(\beta_i) = Z_i$ for $i \in \{K+1, \ldots, K+T\}$ where are $Z_i$'s are chosen randomly from $\bV$.
        \begin{equation}
            u(z) \triangleq \sum_{j\in[K]}X_j\cdot \prod_{k\in [K+T]\setminus\{j\}}\frac{z-\beta_k}{\beta_j-\beta_k} + \sum_{j=K+1}^{K+T}Z_j\cdot \prod_{k\in [K+T]\setminus\{j\}}\frac{z-\beta_k}{\beta_j-\beta_k}
        \end{equation}
        \pause\item Now, select $N$ distinct elements $\alpha_1,\ldots,\alpha_N$ from $\bF$ such that $\{\alpha_i\}_{i=1}^{N} \cap \{\beta_j\}_{j=1}^{N} = \varnothing$ and  encode the input variables by letting $\tilde{X} = u(\alpha_i)$ for $i \in [N]$
        \begin{equation}
            \tilde{X}_i = g_i(X) = u(\alpha_i) = (X_1,\ldots,X_K,Z_{K+1},\ldots,Z_{K+T})\cdot U_i \quad (U \in \bF_q^{(K+T)\times N)})
        \end{equation}
    Where $ U_{ij} \triangleq \prod_{l \in [K+T]\setminus\{i\}}\frac{\alpha_j-\beta_l}{\beta_i-\beta_l}$
    \end{itemize}
\end{frame}
\begin{frame}{Lagrange Coded Computing}{Decoding}
    \begin{itemize}
        \pause\item Each worker $i$ computes $\tilde{Y}_i = f(\tilde{X}_i) = f(u(\alpha_i)$ and sends $\tilde{Y}_i$ to the master.
        \pause\item This composition $f(u(z))$ is also a polynomial with degree $\leq (K+T-1)\deg f$. 
        \pause\item The master obtains $N-S$ evaluations of $f(u(z))$, at most $A$ of which are incorrect.
        \pause\item The master can obtain all coefficients of~$f(u(z))$ by applying Reed-Solomon decoding as $N \geq (K + T - 1) \deg(f) + S + 2A + 1$
        \pause\item Then, the master evaluates it at~$\beta_i$ for every~$i\in[K]$ to obtain~$f(u(\beta_i))=f(X_i)$, %\blue
    \end{itemize}
    Hence, we have shown that the above scheme is $S$-resilient and $A$-secure.
\end{frame}
\begin{frame}{Recent Works and Open Problems}
    \begin{itemize}
        \pause \item Developing  efficient and straggler resilient ``master-less'' systems.
\begin{itemize}
    \pause \item Current state of the art assumes availability of master
    \pause \item All nodes are identical, and no single node may be able to store all the data or perform all the encoding/decoding.
\end{itemize}
\pause  \item Beyond polynomial computations.
\begin{itemize}
    \pause \item Going beyond polynomial computations is a very important and challenging research 
    \pause \item Impacts application domains (eg. machine learning with non-linear threshold functions)
\end{itemize}
\pause \item Application to blockchain systems
\begin{itemize}
    \pause \item Today’s blockchain designs suffer from a trilemma claiming that no blockchain system can simultaneously achieve decentralization, security, and performance scalability.
    \pause \item Coded computing can provide an effective approach for overcoming such barriers.
\end{itemize}
    \end{itemize}
\end{frame}
% \begin{frame}{Blocks}
% \begin{block}{Block Title}
% You can also highlight sections of your presentation in a block, with it's own title
% \end{block}
% \begin{theorem}
% There are separate environments for theorems, examples, definitions and proofs.
% \end{theorem}
% \begin{example}
% Here is an example of an example block.
% \end{example}
% \begin{proof}
% Here is an example of an proof block.
% \end{proof}
% \end{frame}

% % Placing a * after \section means it will not show in the
% % outline or table of contents.
% \section*{Summary}


% \begin{frame}{Summary}
%   \begin{itemize}
%   \pause\item
%     The \alert{first main message} of your talk in one or two lines.
%   \pause\item
%     The \alert{second main message} of your talk in one or two lines.
%   \pause\item
%     Perhaps a \alert{third message}, but not more than that.
%   \end{itemize}
  
%   \begin{itemize}
%   \pause\item
%     Outlook
%     \begin{itemize}
%     \pause\item
%       Something you haven't solved.
%     \pause\item
%       Something else you haven't solved.
%     \end{itemize}
%   \end{itemize}
% \end{frame}



% All of the following is optional and typically not needed. 
\appendix
\section<presentation>*{\appendixname}
\subsection<presentation>*{For Further Reading}
% \begin{frame}[allowframebreaks]
\begin{frame}
  \frametitle<presentation>{For Further Reading}
%   \bibliographystyle{plainurl}    
\bibliographystyle{ieeetr}
  \nocite{*}
  \bibliography{references}
%   \begin{thebibliography}{10}
    
%   \beamertemplatebookbibitems
%   % Start with overview books.

%   \bibitem{Author1990}
%     A.~Author.
%     \newblock {\em Handbook of Everything}.
%     \newblock Some Press, 1990.
 
    
%   \beamertemplatearticlebibitems
%   % Followed by interesting articles. Keep the list short. 

%   \bibitem{Someone2000}
%     S.~Someone.
%     \newblock On this and that.
%     \newblock {\em Journal of This and That}, 2(1):50--100,
%     2000.
%   \end{thebibliography}
\end{frame}

\end{document}


