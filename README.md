# SurvBin Package

The \texttt{SurvBin} package contains three key functions: 
\texttt{lstats} to compute the standardized $\mathcal{L}$-statistic, {\small$\mathrm{\textbf{U}}_{n}^{\omega}(\hat{Q})\big/\sqrt{\widehat{\mathrm{Var}}(\mathrm{\textbf{U}}_{n}^{\omega}(\hat{Q}))}$}; and 
\texttt{bintest} and 
\texttt{survtest}  for the binary and survival statistics, {\small$U_{b,n}\left(\tau_b\right)\big/\hat{\sigma}_b$} and {\small$U_{s,n}(\tau_0,\tau;\hat{Q})\big/\hat{\sigma}_s$}, respectively. 
The \texttt{SurvBin} package also provides the function \texttt{survbinCov} that can be used to calculate the covariance $\sigma_{bs}$; and the functions \texttt{simsurvbin} and \texttt{simsurv}  for simulating bivariate binary and survival and univariate survival data, respectively.