# 统计学的基本知识

````{prf:definition}
:label: def-1.1

一个统计问题所研究的对象的全体称为总体(Population)，而样本(Sample) 是从总体中抽取的一部分个体。在数理统计学中，总体可以用一个随机变量X 及其概率空间$\{\Omega,\mathcal{F}_0,P\}$ 来描述。  

````

· **定义1.2**

一个问题的统计模型(Statistical Model)是指研究该问题时所抽样本的分布。具体而言：  

- 随机向量 $\mathbf{X} = ( X_1, \cdots , X_n)$（样本）: 观测数据 $X_1, \cdots , X_n\sim X.$  
- $\mathbf{X}$ 的概率分布 $( \mathcal{X} , \mathcal{F} , P)$；样本空间 $\mathcal{X}$（sample space）即样本 $\mathbf{X}$ 可能取值的全体；$\mathcal{F}$为样本所发生的所有可测事件组成的$\sigma$代数，P是在可测空间$(\mathcal{X},\mathcal{F})$上定义的一个概率分布。  
- 统计模型$\mathcal{P}=\{P:P$是$(\mathcal{X},\mathcal{F})$上的一个概率分布族$\}$。  
  
· **定义1.3**

对于统计模型$\mathcal{P}$，假如分布族仅依赖于某个参数（或参数向量）$\theta$，即  
$$
\mathcal{P}=\{P_\theta:\theta\in\Theta\}  
$$
其中$\Theta$为参数空间，则称此模型为参数（化）统计模型。  
  
**注**：如果统计模型$\mathcal{P}$中不含参数，则称之为非参数统计模型；半参数统计模型即既含有非参数部分又含有参数部分的统计模型。
