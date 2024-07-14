# 第一章

这是第一章的介绍

In this Section, we state the asymptotic normality of the estimator. First, we need some preliminaries. Let $\rho>0$ and $\mathcal{B}_\rho\left(\boldsymbol{\theta}_0\right)=\left\{\boldsymbol{\theta} \in \Theta \mid\left\|\boldsymbol{\theta}-\boldsymbol{\theta}_0\right\| \leq \rho\right\}$ be a ball around $\boldsymbol{\theta}_0$. Since $\boldsymbol{\theta}_0 \in \Theta$, for sufficiently small $\rho>0, \mathcal{B}_\rho\left(\boldsymbol{\theta}_0\right) \in \Theta$. Let $\mathcal{L}_N$ be one of the two objective functions (22) or (23). For $\hat{\boldsymbol{\theta}}_N \in \mathcal{B}_\rho\left(\boldsymbol{\theta}_0\right)$, the mean value theorem yields:
$$
\left(\int_0^1 \mathbf{H}_{\mathcal{L}_N}\left(\boldsymbol{\theta}_0+t\left(\hat{\boldsymbol{\theta}}_N-\boldsymbol{\theta}_0\right)\right) \mathrm{d} t\right)\left(\hat{\boldsymbol{\theta}}_N-\boldsymbol{\theta}_0\right)=-\nabla \mathcal{L}_N\left(\boldsymbol{\theta}_0\right) .
$$

With $\boldsymbol{\varsigma}:=\operatorname{vech}\left(\boldsymbol{\Sigma} \boldsymbol{\Sigma}^{\top}\right)=\left(\left[\boldsymbol{\Sigma} \boldsymbol{\Sigma}^{\top}\right]_{11},\left[\boldsymbol{\Sigma} \boldsymbol{\Sigma}^{\top}\right]_{12},\left[\boldsymbol{\Sigma} \boldsymbol{\Sigma}^{\top}\right]_{22}, \ldots,\left[\boldsymbol{\Sigma} \boldsymbol{\Sigma}^{\top}\right]_{1 d}, \ldots,\left[\boldsymbol{\Sigma} \boldsymbol{\Sigma}^{\top}\right]_{d d}\right)$, we half-vectorize $\boldsymbol{\Sigma} \boldsymbol{\Sigma}^{\top}$ to avoid working with tensors when computing derivatives with respect to $\boldsymbol{\Sigma} \boldsymbol{\Sigma}^{\top}$. Since $\boldsymbol{\Sigma} \boldsymbol{\Sigma}^{\top}$ is a symmetric $d \times d$ matrix, $\boldsymbol{\varsigma}$ is of dimension $s=d(d+1) / 2$. For a diagonal matrix, instead of a half-vectorization, we use $\boldsymbol{\varsigma}:=\operatorname{diag}\left(\boldsymbol{\Sigma} \boldsymbol{\Sigma}^{\top}\right)$. Define:
$$
\begin{gathered}
\mathbf{C}_N(\boldsymbol{\theta}):=\left[\begin{array}{cc}
\frac{1}{N h} \partial_{\boldsymbol{\beta} \boldsymbol{\beta}} \mathcal{L}_N(\boldsymbol{\theta}) & \frac{1}{N \sqrt{h}} \partial_{\boldsymbol{\beta} \boldsymbol{\varsigma}} \mathcal{L}_N(\boldsymbol{\theta}) \\
\frac{1}{N \sqrt{h}} \partial_{\boldsymbol{\beta} \boldsymbol{\varsigma}} \mathcal{L}_N(\boldsymbol{\theta}) & \frac{1}{N} \partial_{\boldsymbol{\varsigma} \boldsymbol{\varsigma}} \mathcal{L}_N(\boldsymbol{\theta})
\end{array}\right], \\
\mathbf{s}_N:=\left[\begin{array}{c}
\sqrt{N h}\left(\hat{\boldsymbol{\beta}}_N-\boldsymbol{\beta}_0\right) \\
\sqrt{N}\left(\hat{\boldsymbol{\varsigma}}_N-\boldsymbol{\varsigma}_0\right)
\end{array}\right], \quad \boldsymbol{\lambda}_N:=\left[\begin{array}{c}
-\frac{1}{\sqrt{N h}} \partial_{\boldsymbol{\beta}} \mathcal{L}_N\left(\boldsymbol{\theta}_0\right) \\
-\frac{1}{\sqrt{N}} \partial_{\boldsymbol{\varsigma}} \mathcal{L}_N\left(\boldsymbol{\theta}_0\right)
\end{array}\right],
\end{gathered}
$$
and $\mathbf{D}_N:=\int_0^1 \mathbf{C}_N\left(\boldsymbol{\theta}_0+t\left(\hat{\boldsymbol{\theta}}_N-\boldsymbol{\theta}_0\right)\right) \mathrm{d} t$. Then, (24) is equivalent to $\mathbf{D}_N \mathbf{s}_N=\boldsymbol{\lambda}_N$. Let:
$$
\mathbf{C}\left(\boldsymbol{\theta}_0\right):=\left[\begin{array}{cc}
\mathbf{C}_\beta\left(\boldsymbol{\theta}_0\right) & \mathbf{0}_{r \times s} \\
\mathbf{0}_{s \times r} & \mathbf{C}_{\varsigma}\left(\boldsymbol{\theta}_0\right)
\end{array}\right]
$$
where:
$$
\begin{aligned}
{\left[\mathbf{C}_\beta\left(\boldsymbol{\theta}_0\right)\right]_{i_1, i_2} } & :=\int\left(\partial_{\beta_{i_1}} \mathbf{F}_0(\mathbf{x})\right)^{\top}\left(\boldsymbol{\Sigma} \boldsymbol{\Sigma}_0^{\top}\right)^{-1}\left(\partial_{\beta_{i_2}} \mathbf{F}_0(\mathbf{x})\right) \mathrm{d} \nu_0(\mathbf{x}), 1 \leq i_1, i_2 \leq r \\
{\left[\mathbf{C}_{\varsigma}\left(\boldsymbol{\theta}_0\right)\right]_{j_1, j_2} } & :=\frac{1}{2} \operatorname{Tr}\left(\left(\partial \varsigma_{j_1} \boldsymbol{\Sigma} \boldsymbol{\Sigma}_0^{\top}\right)\left(\boldsymbol{\Sigma} \boldsymbol{\Sigma}_0^{\top}\right)^{-1}\left(\partial \varsigma_{j_2} \boldsymbol{\Sigma} \boldsymbol{\Sigma}_0^{\top}\right)\left(\boldsymbol{\Sigma} \boldsymbol{\Sigma}_0^{\top}\right)^{-1}\right), 1 \leq j_1, j_2 \leq s
\end{aligned}
$$

Now, we state the theorem for asymptotic normality, whose proof is in Section 7.4.
Theorem 5.2 Let Assumptions (A1)-(A6) hold, $\mathbf{X}$ be the solution of (1), and $\widehat{\boldsymbol{\theta}}_N=\left(\widehat{\boldsymbol{\beta}}_N, \widehat{\boldsymbol{\varsigma}}_N\right)$ be the estimator that minimizes one of the objective functions (22) or (23). If $\boldsymbol{\theta}_0 \in \Theta, \mathbf{C}\left(\boldsymbol{\theta}_0\right)$ is positive definite, $h \rightarrow 0, N h \rightarrow \infty$, and $N h^2 \rightarrow 0$, then,
$$
\left[\begin{array}{c}
\sqrt{N h}\left(\hat{\boldsymbol{\beta}}_N-\boldsymbol{\beta}_0\right) \\
\sqrt{N}\left(\hat{\boldsymbol{\varsigma}}_N-\boldsymbol{\varsigma}_0\right)
\end{array}\right] \xrightarrow{d} \mathcal{N}\left(\mathbf{0}, \mathbf{C}^{-1}\left(\boldsymbol{\theta}_0\right)\right),
$$
under $\mathbb{P}_{\boldsymbol{\theta}_0}$.
The estimator of the diffusion parameter converges faster than the estimator of the drift parameter. Gobet (2002) showed that for a discretely sampled SDE model, the optimal convergence rates for the drift and diffusion parameters are $1 / \sqrt{N h}$ and $1 / \sqrt{N}$, respectively. Thus, our estimators reach optimal rates. Moreover, the estimators are asymptotically efficient since $\mathbf{C}$ is the Fisher information matrix for the corresponding continuous-time diffusion (see Kessler (1997), Gobet (2002)). Finally, since the asymptotic correlation is zero between the drift and diffusion estimators , they are asymptotically independent.
