<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build
Status](https://api.travis-ci.org/kassambara/rstatix.png)](https://travis-ci.org/kassambara/rstatix)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/rstatix)](https://cran.r-project.org/package=rstatix)
[![CRAN
Checks](https://cranchecks.info/badges/summary/rstatix)](https://cran.r-project.org/web/checks/check_results_rstatix.html)
[![Downloads](https://cranlogs.r-pkg.org/badges/rstatix)](https://cran.r-project.org/package=rstatix)
[![Total
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/rstatix?color=orange)](https://cran.r-project.org/package=rstatix)

rstatix
=======

Provides a simple and intuitive pipe-friendly framework, coherent with
the ‘tidyverse’ design philosophy, for performing basic statistical
tests, including t-test, Wilcoxon test, ANOVA, Kruskal-Wallis and
correlation analyses.

The output of each test is automatically transformed into a tidy data
frame to facilitate visualization.

Additional functions are available for reshaping, reordering,
manipulating and visualizing correlation matrix. Functions are also
included to facilitate the analysis of factorial experiments, including
purely ‘within-Ss’ designs (repeated measures), purely ‘between-Ss’
designs, and mixed ‘within-and-between-Ss’ designs.

It’s also possible to compute several effect size metrics, including
“eta squared” for ANOVA, “Cohen’s d” for t-test and “Cramer’s V” for the
association between categorical variables. The package contains helper
functions for identifying univariate and multivariate outliers,
assessing normality and homogeneity of variances.

核心函数
-------------

### Descriptive statistics

-   `get_summary_stats()`: 计算一个或多个数值变量的摘要统计信息。可以处理分组数据.
-   `freq_table()`: 计算分类变量频率表.
-   `get_mode()`: 计算向量的模式，这是最常见的值.
-   `identify_outliers()`: 用箱线图方法检测单变量异常值.
-   `mahalanobis_distance()`: 计算马氏距离和标志多元离群值.
-   `shapiro_test()` and `mshapiro_test()`: 单变量和多变量Shapiro-Wilk正态性检验.

### 均值比较

-   `t_test()`: 进行单样本、两样本和成对t检验
-   `wilcox_test()`: 进行单样本、两样本和成对的Wilcoxon测试
-   `sign_test()`: 进行符号测试, 以确定成对或匹配的观察值之间是否存在中位数差异.
-   `anova_test()`: 使用`car::Anova()`进行不同类型方差分析的易用包装, 包括__单次测量方差分析__, __重复测量方差分析__以及__混合方差分析__.
-   `get_anova_test_table()`: 从`anova_test()`中提取方差表。 可在受试者内（重复测量）设计中自动应用球度校正.
    `- welch_anova_test()`: 韦尔奇 (Welch) 单因素方差分析. A pipe-friendly
    wrapper around the base function `stats::oneway.test()`. 在方差齐性假设被违反的情况下，这是标准单因素方差分析的一种替代方法.
-   `kruskal_test()`: 进行kruskal-wallis秩和检验
-   `friedman_test()`: 提供了一个管道友好的框架来执行Friedman秩和检验，它是单向重复测量方差分析的非参数替代方法.
-   `get_comparisons()`: 创建组间可能的成对比较列表.
-   `get_pvalue_position`: 使用ggplot2自动计算绘制显著性的p值位置.

### 增强 R 的方差分析计算

-   `factorial_design()`: build factorial design for easily computing
    ANOVA using the `car::Anova()` function. This might be very useful
    for repeated measures ANOVA, which is hard to set up with the `car`
    package.
-   `anova_summary()`: Create beautiful summary tables of ANOVA test
    results obtained from either `car::Anova()` or `stats::aov()`. 结果包括方差分析表、广义效应大小和一些假设检验，如重复测量方差分析中的Mauchly球度检验.

### 事后 (post-hoc) 分析

-   `tukey_hsd()`: 执行tukey事后测试。可以处理不同的输入格式：aov，lm，formula.
-   `dunn_test()`: 根据Kruskal-Wallis检验计算多对比较.
-   `games_howell_test()`: 进行博弈豪厄尔检验，当方差齐性假设被违背时，用来比较所有可能的群体差异组合.
-   `emmeans_test()`: pipe-friendly wrapper arround `emmeans` function
    to perform pairwise comparisons of estimated marginal means. 用于ANOVA/ANCOVA测试后的事后分析.

### 比较比例数据

-   `prop_test()`, `pairwise_prop_test()` 和 `row_wise_prop_test()`.
    执行一个样本和两个样本的比例z检验. Wrappers
    around the R base function `prop.test()` but have the advantage of
    performing pairwise and row-wise z-test of two proportions, the
    post-hoc tests following a significant chi-square test of
    homogeneity for 2xc and rx2 contingency tables.
-   `fisher_test()`, `pairwise_fisher_test()` and
    `row_wise_fisher_test()`: 计数数据的Fisher精确检验.
    Wrappers around the R base function `fisher.test()` but have the
    advantage of performing pairwise and row-wise fisher tests, the
    post-hoc tests following a significant chi-square test of
    homogeneity for 2xc and rx2 contingency tables.
-   `chisq_test()`, `pairwise_chisq_gof_test()`,
    `pairwise_chisq_test_against_p()`: 进行卡方检验，包括拟合优度、同质性和独立性检验.
-   `binom_test()`, `pairwise_binom_test()`,
    `pairwise_binom_test_against_p()`: 执行精确的二项式检验和配对比较，然后进行显著的精确多项式检验.
    替代卡方检验的拟合优度检验当样本.
-   `multinom_test()`: 进行精确的多项式检验. 当样本量较小时，可替代卡方检验的拟合优度检验.
-   `mcnemar_test()`: 执行McNemar卡方检验来比较配对比例. 提供多个组之间的成对比较.
-   `cochran_qtest()`: McNemar卡方检验用于比较两个以上配对比例的扩展.
-   `prop_trend_test()`: 按比例对趋势进行卡方检验。这种测试也被称为Cochran-Armitage趋势测试.

### 比较方差

-   `levene_test()`: 管道友好的框架，可以轻松地计算各组方差齐性的Levene检验。处理分组数据.
-   `box_m()`: 协方差矩阵齐性的Box M检验

### 效应量

-   `cohens_d()`: t检验中cohen的效应大小度量.
-   `wilcox_effsize()`: Compute Wilcoxon effect size (r).
-   `eta_squared()` and `partial_eta_squared()`: Compute effect size for
    ANOVA.
-   `kruskal_effsize()`: Compute the effect size for Kruskal-Wallis test
    as the eta squared based on the H-statistic.
-   `friedman_effsize()`: Compute the effect size of Friedman test using
    the Kendall’s W value.
-   `cramer_v()`: Compute Cramer’s V, which measures the strength of the
    association between categorical variables.

### 相关性分析

**Computing correlation**:

-   `cor_test()`: correlation test between two or more variables using
    Pearson, Spearman or Kendall methods.
-   `cor_mat()`: compute correlation matrix with p-values. Returns a
    data frame containing the matrix of the correlation coefficients.
    The output has an attribute named “pvalue”, which contains the
    matrix of the correlation test p-values.
-   `cor_get_pval()`: extract a correlation matrix p-values from an
    object of class `cor_mat()`.
-   `cor_pmat()`: compute the correlation matrix, but returns only the
    p-values of the correlation tests.
-   `as_cor_mat()`: convert a `cor_test` object into a correlation
    matrix format.

**Reshaping correlation matrix**:

-   `cor_reorder()`: reorder correlation matrix, according to the
    coefficients, using the hierarchical clustering method.
-   `cor_gather()`: takes a correlation matrix and collapses (or melt)
    it into long format data frame (paired list)
-   `cor_spread()`: spread a long correlation data frame into wide
    format (correlation matrix).

**Subsetting correlation matrix**:

-   `cor_select()`: subset a correlation matrix by selecting variables
    of interest.
-   `pull_triangle()`, `pull_upper_triangle()`, `pull_lower_triangle()`:
    pull upper and lower triangular parts of a (correlation) matrix.
-   `replace_triangle()`, `replace_upper_triangle()`,
    `replace_lower_triangle()`: replace upper and lower triangular parts
    of a (correlation) matrix.

**Visualizing correlation matrix**:

-   `cor_as_symbols()`: replaces the correlation coefficients, in a
    matrix, by symbols according to the value.
-   `cor_plot()`: visualize correlation matrix using base plot.
-   `cor_mark_significant()`: add significance levels to a correlation
    matrix.

### Adjusting p-values, formatting and adding significance symbols

-   `adjust_pvalue()`: add an adjusted p-values column to a data frame
    containing statistical test p-values
-   `add_significance()`: add a column containing the p-value
    significance level
-   `p_round(), p_format(), p_mark_significant()`: rounding and
    formatting p-values

### Extract information from statistical tests

Extract information from statistical test results. Useful for labelling
plots with test outputs.

-   `get_pwc_label()`: Extract label from pairwise comparisons.
-   `get_test_label()`: Extract label from statistical tests.
-   `create_test_label()`: Create labels from user specified test
    results.

### Data manipulation helper functions

These functions are internally used in the `rstatix` and in the `ggpubr`
R package to make it easy to program with tidyverse packages using non
standard evaluation.

-   `df_select()`, `df_arrange()`, `df_group_by()`: wrappers arround
    dplyr functions for supporting standard and non standard
    evaluations.
-   `df_nest_by()`: Nest a tibble data frame using grouping
    specification. Supports standard and non standard evaluations.
-   `df_split_by()`: Split a data frame by groups into subsets or data
    panel. Very similar to the function `df_nest_by()`. The only
    difference is that, it adds labels to each data subset. Labels are
    the combination of the grouping variable levels.
-   `df_unite()`: Unite multiple columns into one.
-   `df_unite_factors()`: Unite factor columns. First, order factors
    levels then merge them into one column. The output column is a
    factor.
-   `df_label_both()`, `df_label_value()`: functions to label data
    frames rows by by one or multiple grouping variables.
-   `df_get_var_names()`: Returns user specified variable names.
    Supports standard and non standard evaluation.

### Others

-   `doo()`: alternative to dplyr::do for doing anything. Technically it
    uses `nest() + mutate() + map()` to apply arbitrary computation to a
    grouped data frame.
-   `sample_n_by()`: sample n rows by group from a table
-   `convert_as_factor(), set_ref_level(), reorder_levels()`: Provides
    pipe-friendly functions to convert simultaneously multiple variables
    into a factor variable.
-   `make_clean_names()`: Pipe-friendly function to make syntactically
    valid column names (for input data frame) or names (for input
    vector).
-   `counts_to_cases()`: converts a contingency table or a data frame of
    counts into a data frame of individual observations.

包的安装和加载
------------------------

-   Install the latest developmental version from
    [GitHub](https://github.com/kassambara/rstatix) as follow:

``` r
if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/rstatix")
```

-   Or install from [CRAN](https://cran.r-project.org/package=ggpubr) as
    follow:

``` r
install.packages("rstatix")
```

-   Loading packages

``` r
library(rstatix)  
library(ggpubr)  # For easy data-visualization
```

描述性统计
----------------------

``` r
# Summary statistics of some selected variables
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
iris %>% 
  get_summary_stats(Sepal.Length, Sepal.Width, type = "common")
#> # A tibble: 2 x 10
#>   variable         n   min   max median   iqr  mean    sd    se    ci
#>   <chr>        <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1 Sepal.Length   150   4.3   7.9    5.8   1.3  5.84 0.828 0.068 0.134
#> 2 Sepal.Width    150   2     4.4    3     0.5  3.06 0.436 0.036 0.07

# Whole data frame
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
iris %>% get_summary_stats(type = "common")
#> # A tibble: 4 x 10
#>   variable         n   min   max median   iqr  mean    sd    se    ci
#>   <chr>        <dbl> <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
#> 1 Petal.Length   150   1     6.9   4.35   3.5  3.76 1.76  0.144 0.285
#> 2 Petal.Width    150   0.1   2.5   1.3    1.5  1.20 0.762 0.062 0.123
#> 3 Sepal.Length   150   4.3   7.9   5.8    1.3  5.84 0.828 0.068 0.134
#> 4 Sepal.Width    150   2     4.4   3      0.5  3.06 0.436 0.036 0.07


# Grouped data
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
iris %>%
  group_by(Species) %>% 
  get_summary_stats(Sepal.Length, type = "mean_sd")
#> # A tibble: 3 x 5
#>   Species    variable         n  mean    sd
#>   <fct>      <chr>        <dbl> <dbl> <dbl>
#> 1 setosa     Sepal.Length    50  5.01 0.352
#> 2 versicolor Sepal.Length    50  5.94 0.516
#> 3 virginica  Sepal.Length    50  6.59 0.636
```

两均值比较
-------------------

要比较两组的平均值，可以使用函数`t_test()` (参数) 或`wilcox_test()` (非参). 在下面的例子中，将说明t检验.

### 数据

Preparing the demo data set:

``` r
df <- ToothGrowth
df$dose <- as.factor(df$dose)
head(df)
#>    len supp dose
#> 1  4.2   VC  0.5
#> 2 11.5   VC  0.5
#> 3  7.3   VC  0.5
#> 4  5.8   VC  0.5
#> 5  6.4   VC  0.5
#> 6 10.0   VC  0.5
```

### 单样本检验

The one-sample test is used to compare the mean of one sample to a known
standard (or theoretical / hypothetical) mean (`mu`).

``` r
df %>% t_test(len ~ 1, mu = 0)
#> # A tibble: 1 x 7
#>   .y.   group1 group2         n statistic    df        p
#> * <chr> <chr>  <chr>      <int>     <dbl> <dbl>    <dbl>
#> 1 len   1      null model    60      19.1    59 6.94e-27
# One-sample test of each dose level
df %>% 
  group_by(dose) %>%
  t_test(len ~ 1, mu = 0)
#> # A tibble: 3 x 8
#>   dose  .y.   group1 group2         n statistic    df        p
#> * <fct> <chr> <chr>  <chr>      <int>     <dbl> <dbl>    <dbl>
#> 1 0.5   len   1      null model    20      10.5    19 2.24e- 9
#> 2 1     len   1      null model    20      20.0    19 3.22e-14
#> 3 2     len   1      null model    20      30.9    19 1.03e-17
```

### 两组独立样本检验

-   创建带有p值的箱型图

``` r
# T-test
stat.test <- df %>% 
  t_test(len ~ supp, paired = FALSE) 
stat.test
#> # A tibble: 1 x 8
#>   .y.   group1 group2    n1    n2 statistic    df      p
#> * <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl>  <dbl>
#> 1 len   OJ     VC        30    30      1.92  55.3 0.0606

# Create a box plot
p <- ggboxplot(
  df, x = "supp", y = "len", 
  color = "supp", palette = "jco", ylim = c(0,40)
  )
# Add the p-value manually
p + stat_pvalue_manual(stat.test, label = "p", y.position = 35)
```

![](tools/README-unpaired-two-sample-t-test-1.png)

-   Customize labels using [glue
    expression](https://github.com/tidyverse/glue):

``` r
p +stat_pvalue_manual(stat.test, label = "T-test, p = {p}", 
                      y.position = 36)
```

![](tools/README-custoize-p-value-labels-1.png)

-   数据分组: compare supp levels after grouping the data by “dose”

``` r
# Statistical test
stat.test <- df %>%
  group_by(dose) %>%
  t_test(len ~ supp) %>%
  adjust_pvalue() %>%
  add_significance("p.adj")
stat.test
#> # A tibble: 3 x 11
#>   dose  .y.   group1 group2    n1    n2 statistic    df       p   p.adj
#>   <fct> <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl>   <dbl>   <dbl>
#> 1 0.5   len   OJ     VC        10    10    3.17    15.0 0.00636 0.0127 
#> 2 1     len   OJ     VC        10    10    4.03    15.4 0.00104 0.00312
#> 3 2     len   OJ     VC        10    10   -0.0461  14.0 0.964   0.964  
#> # … with 1 more variable: p.adj.signif <chr>

# Visualization
ggboxplot(
  df, x = "supp", y = "len",
  color = "supp", palette = "jco", facet.by = "dose",
  ylim = c(0, 40)
  ) +
  stat_pvalue_manual(stat.test, label = "p.adj", y.position = 35)
```

![](tools/README-grouped-two-sample-t-test-1.png)

### 成对样本检验

``` r
# T-test
stat.test <- df %>% 
  t_test(len ~ supp, paired = TRUE) 
stat.test
#> # A tibble: 1 x 8
#>   .y.   group1 group2    n1    n2 statistic    df       p
#> * <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl>   <dbl>
#> 1 len   OJ     VC        30    30      3.30    29 0.00255

# Box plot
p <- ggpaired(
  df, x = "supp", y = "len", color = "supp", palette = "jco", 
  line.color = "gray", line.size = 0.4, ylim = c(0, 40)
  )
p + stat_pvalue_manual(stat.test, label = "p", y.position = 36)
```

![](tools/README-paired-t-test-1.png)

### 多重比较

-   Pairwise comparisons: if the grouping variable contains more than
    two categories, a pairwise comparison is automatically performed.

``` r
# Pairwise t-test
pairwise.test <- df %>% t_test(len ~ dose)
pairwise.test
#> # A tibble: 3 x 10
#>   .y.   group1 group2    n1    n2 statistic    df        p    p.adj p.adj.signif
#> * <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl>    <dbl>    <dbl> <chr>       
#> 1 len   0.5    1         20    20     -6.48  38.0 1.27e- 7 2.54e- 7 ****        
#> 2 len   0.5    2         20    20    -11.8   36.9 4.40e-14 1.32e-13 ****        
#> 3 len   1      2         20    20     -4.90  37.1 1.91e- 5 1.91e- 5 ****
# Box plot
ggboxplot(df, x = "dose", y = "len")+
  stat_pvalue_manual(
    pairwise.test, label = "p.adj", 
    y.position = c(29, 35, 39)
    )
```

![](tools/README-pairwise-comparisons-1.png)

-   Multiple pairwise comparisons against reference group: each level is
    compared to the ref group

``` r
# Comparison against reference group
#::::::::::::::::::::::::::::::::::::::::
# T-test: each level is compared to the ref group
stat.test <- df %>% t_test(len ~ dose, ref.group = "0.5")
stat.test
#> # A tibble: 2 x 10
#>   .y.   group1 group2    n1    n2 statistic    df        p    p.adj p.adj.signif
#> * <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl>    <dbl>    <dbl> <chr>       
#> 1 len   0.5    1         20    20     -6.48  38.0 1.27e- 7 1.27e- 7 ****        
#> 2 len   0.5    2         20    20    -11.8   36.9 4.40e-14 8.80e-14 ****
# Box plot
ggboxplot(df, x = "dose", y = "len", ylim = c(0, 40)) +
  stat_pvalue_manual(
    stat.test, label = "p.adj.signif", 
    y.position = c(29, 35)
    )
```

![](tools/README-comaprison-against-reference-group-1.png)

``` r
# Remove bracket
ggboxplot(df, x = "dose", y = "len", ylim = c(0, 40)) +
  stat_pvalue_manual(
    stat.test, label = "p.adj.signif", 
    y.position = c(29, 35),
    remove.bracket = TRUE
    )
```

![](tools/README-comaprison-against-reference-group-2.png)

-   Multiple pairwise comparisons against all (base-mean): Comparison of
    each group against base-mean.

``` r
# T-test
stat.test <- df %>% t_test(len ~ dose, ref.group = "all")
stat.test
#> # A tibble: 3 x 10
#>   .y.   group1 group2    n1    n2 statistic    df         p   p.adj p.adj.signif
#> * <chr> <chr>  <chr>  <int> <int>     <dbl> <dbl>     <dbl>   <dbl> <chr>       
#> 1 len   all    0.5       60    20     5.82   56.4   2.90e-7 8.70e-7 ****        
#> 2 len   all    1         60    20    -0.660  57.5   5.12e-1 5.12e-1 ns          
#> 3 len   all    2         60    20    -5.61   66.5   4.25e-7 8.70e-7 ****
# Box plot with horizontal mean line
ggboxplot(df, x = "dose", y = "len") +
  stat_pvalue_manual(
    stat.test, label = "p.adj.signif", 
    y.position = 35,
    remove.bracket = TRUE
    ) +
  geom_hline(yintercept = mean(df$len), linetype = 2)
```

![](tools/README-comparison-against-base-mean-1.png)

方差分析
----------

``` r
# 单因素方差分析
#:::::::::::::::::::::::::::::::::::::::::
df %>% anova_test(len ~ dose)
#> ANOVA Table (type II tests)
#> 
#>   Effect DFn DFd      F        p p<.05   ges
#> 1   dose   2  57 67.416 9.53e-16     * 0.703

# 双因素方差分析
#:::::::::::::::::::::::::::::::::::::::::
df %>% anova_test(len ~ supp*dose)
#> ANOVA Table (type II tests)
#> 
#>      Effect DFn DFd      F        p p<.05   ges
#> 1      supp   1  54 15.572 2.31e-04     * 0.224
#> 2      dose   2  54 92.000 4.05e-18     * 0.773
#> 3 supp:dose   2  54  4.107 2.20e-02     * 0.132

# 双因素重复测量方差分析
#:::::::::::::::::::::::::::::::::::::::::
df$id <- rep(1:10, 6) # Add individuals id
# Use formula
# df %>% anova_test(len ~ supp*dose + Error(id/(supp*dose)))
# or use character vector
df %>% anova_test(dv = len, wid = id, within = c(supp, dose))
#> ANOVA Table (type III tests)
#> 
#> $ANOVA
#>      Effect DFn DFd       F        p p<.05   ges
#> 1      supp   1   9  34.866 2.28e-04     * 0.224
#> 2      dose   2  18 106.470 1.06e-10     * 0.773
#> 3 supp:dose   2  18   2.534 1.07e-01       0.132
#> 
#> $`Mauchly's Test for Sphericity`
#>      Effect     W     p p<.05
#> 1      dose 0.807 0.425      
#> 2 supp:dose 0.934 0.761      
#> 
#> $`Sphericity Corrections`
#>      Effect   GGe      DF[GG]    p[GG] p[GG]<.05   HFe      DF[HF]    p[HF]
#> 1      dose 0.838 1.68, 15.09 2.79e-09         * 1.008 2.02, 18.15 1.06e-10
#> 2 supp:dose 0.938 1.88, 16.88 1.12e-01           1.176 2.35, 21.17 1.07e-01
#>   p[HF]<.05
#> 1         *
#> 2

# 使用模型作为参数
#:::::::::::::::::::::::::::::::::::::::::
.my.model <- lm(yield ~ block + N*P*K, npk)
anova_test(.my.model)
#> ANOVA Table (type II tests)
#> 
#>   Effect DFn DFd      F     p p<.05   ges
#> 1  block   4  12  4.959 0.014     * 0.623
#> 2      N   1  12 12.259 0.004     * 0.505
#> 3      P   1  12  0.544 0.475       0.043
#> 4      K   1  12  6.166 0.029     * 0.339
#> 5    N:P   1  12  1.378 0.263       0.103
#> 6    N:K   1  12  2.146 0.169       0.152
#> 7    P:K   1  12  0.031 0.863       0.003
#> 8  N:P:K   0  12     NA    NA  <NA>    NA
```

相关性分析
-----------------

``` r
# 数据准备
mydata <- mtcars %>% 
  select(mpg, disp, hp, drat, wt, qsec)
head(mydata, 3)
#>                mpg disp  hp drat    wt  qsec
#> Mazda RX4     21.0  160 110 3.90 2.620 16.46
#> Mazda RX4 Wag 21.0  160 110 3.90 2.875 17.02
#> Datsun 710    22.8  108  93 3.85 2.320 18.61

# 两个因子的相关性分析
mydata %>% cor_test(wt, mpg, method = "pearson")
#> # A tibble: 1 x 8
#>   var1  var2    cor statistic        p conf.low conf.high method 
#>   <chr> <chr> <dbl>     <dbl>    <dbl>    <dbl>     <dbl> <chr>  
#> 1 wt    mpg   -0.87     -9.56 1.29e-10   -0.934    -0.744 Pearson

# Correlation of one variable against all
mydata %>% cor_test(mpg, method = "pearson")
#> # A tibble: 5 x 8
#>   var1  var2    cor statistic        p conf.low conf.high method 
#>   <chr> <chr> <dbl>     <dbl>    <dbl>    <dbl>     <dbl> <chr>  
#> 1 mpg   disp  -0.85     -8.75 9.38e-10  -0.923     -0.708 Pearson
#> 2 mpg   hp    -0.78     -6.74 1.79e- 7  -0.885     -0.586 Pearson
#> 3 mpg   drat   0.68      5.10 1.78e- 5   0.436      0.832 Pearson
#> 4 mpg   wt    -0.87     -9.56 1.29e-10  -0.934     -0.744 Pearson
#> 5 mpg   qsec   0.42      2.53 1.71e- 2   0.0820     0.670 Pearson

# Pairwise correlation test between all variables
mydata %>% cor_test(method = "pearson")
#> # A tibble: 36 x 8
#>    var1  var2    cor statistic        p conf.low conf.high method 
#>    <chr> <chr> <dbl>     <dbl>    <dbl>    <dbl>     <dbl> <chr>  
#>  1 mpg   mpg    1       Inf    0.         1          1     Pearson
#>  2 mpg   disp  -0.85     -8.75 9.38e-10  -0.923     -0.708 Pearson
#>  3 mpg   hp    -0.78     -6.74 1.79e- 7  -0.885     -0.586 Pearson
#>  4 mpg   drat   0.68      5.10 1.78e- 5   0.436      0.832 Pearson
#>  5 mpg   wt    -0.87     -9.56 1.29e-10  -0.934     -0.744 Pearson
#>  6 mpg   qsec   0.42      2.53 1.71e- 2   0.0820     0.670 Pearson
#>  7 disp  mpg   -0.85     -8.75 9.38e-10  -0.923     -0.708 Pearson
#>  8 disp  disp   1       Inf    0.         1          1     Pearson
#>  9 disp  hp     0.79      7.08 7.14e- 8   0.611      0.893 Pearson
#> 10 disp  drat  -0.71     -5.53 5.28e- 6  -0.849     -0.481 Pearson
#> # … with 26 more rows
```

相关性矩阵
------------------

``` r
# Compute correlation matrix
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
cor.mat <- mydata %>% cor_mat()
cor.mat
#> # A tibble: 6 x 7
#>   rowname   mpg  disp    hp   drat    wt   qsec
#> * <chr>   <dbl> <dbl> <dbl>  <dbl> <dbl>  <dbl>
#> 1 mpg      1    -0.85 -0.78  0.68  -0.87  0.42 
#> 2 disp    -0.85  1     0.79 -0.71   0.89 -0.43 
#> 3 hp      -0.78  0.79  1    -0.45   0.66 -0.71 
#> 4 drat     0.68 -0.71 -0.45  1     -0.71  0.091
#> 5 wt      -0.87  0.89  0.66 -0.71   1    -0.17 
#> 6 qsec     0.42 -0.43 -0.71  0.091 -0.17  1

# Show the significance levels
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
cor.mat %>% cor_get_pval()
#> # A tibble: 6 x 7
#>   rowname      mpg     disp           hp       drat        wt       qsec
#>   <chr>      <dbl>    <dbl>        <dbl>      <dbl>     <dbl>      <dbl>
#> 1 mpg     0.       9.38e-10 0.000000179  0.0000178  1.29e- 10 0.0171    
#> 2 disp    9.38e-10 0.       0.0000000714 0.00000528 1.22e- 11 0.0131    
#> 3 hp      1.79e- 7 7.14e- 8 0            0.00999    4.15e-  5 0.00000577
#> 4 drat    1.78e- 5 5.28e- 6 0.00999      0          4.78e-  6 0.62      
#> 5 wt      1.29e-10 1.22e-11 0.0000415    0.00000478 2.27e-236 0.339     
#> 6 qsec    1.71e- 2 1.31e- 2 0.00000577   0.62       3.39e-  1 0

# Replacing correlation coefficients by symbols
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
cor.mat %>%
  cor_as_symbols() %>%
  pull_lower_triangle()
#>   rowname mpg disp hp drat wt qsec
#> 1     mpg                         
#> 2    disp   *                     
#> 3      hp   *    *                
#> 4    drat   +    +  .             
#> 5      wt   *    *  +    +        
#> 6    qsec   .    .  +

# Mark significant correlations
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
cor.mat %>%
  cor_mark_significant()
#>   rowname       mpg      disp        hp      drat    wt qsec
#> 1     mpg                                                   
#> 2    disp -0.85****                                         
#> 3      hp -0.78****  0.79****                               
#> 4    drat  0.68**** -0.71****   -0.45**                     
#> 5      wt -0.87****  0.89****  0.66**** -0.71****           
#> 6    qsec     0.42*    -0.43* -0.71****     0.091 -0.17


# Draw correlogram using R base plot
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
cor.mat %>%
  cor_reorder() %>%
  pull_lower_triangle() %>% 
  cor_plot()
```

![](tools/README-unnamed-chunk-10-1.png)

Related articles
----------------

-   [How to Add P-Values onto Basic
    GGPLOTS](https://www.datanovia.com/en/blog/how-to-add-p-values-onto-basic-ggplots/)
-   [How to Add Adjusted P-values to a Multi-Panel
    GGPlot](https://www.datanovia.com/en/blog/ggpubr-how-to-add-adjusted-p-values-to-a-multi-panel-ggplot/)
-   [How to Add P-values to GGPLOT
    Facets](https://www.datanovia.com/en/blog/how-to-add-p-values-to-ggplot-facets/)
-   [How to Add P-Values Generated Elsewhere to a
    GGPLOT](https://www.datanovia.com/en/blog/ggpubr-how-to-add-p-values-generated-elsewhere-to-a-ggplot/)
-   [How to Add P-Values onto a Grouped GGPLOT using the GGPUBR R
    Package](https://www.datanovia.com/en/blog/how-to-add-p-values-onto-a-grouped-ggplot-using-the-ggpubr-r-package/)
-   [How to Create Stacked Bar Plots with Error Bars and
    P-values](https://www.datanovia.com/en/blog/how-to-create-stacked-bar-plots-with-error-bars-and-p-values/)
-   [How to Add P-Values onto Horizontal
    GGPLOTS](https://www.datanovia.com/en/blog/how-to-add-p-values-onto-horizontal-ggplots/)
-   [Add P-values and Significance Levels to
    ggplots](http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/)
-   [Comparing Means of Two Groups in
    R](https://www.datanovia.com/en/courses/comparing-means-of-two-groups-in-r/)
    -   [T-test in R](https://www.datanovia.com/en/lessons/t-test-in-r/)
    -   [Wilcoxon Test in
        R](https://www.datanovia.com/en/lessons/wilcoxon-test-in-r/)
    -   [Sign Test in
        R](https://www.datanovia.com/en/lessons/sign-test-in-r/)
-   [Comparing Multiple Means in
    R](https://www.datanovia.com/en/courses/comparing-multiple-means-in-r/)
    -   [ANOVA in R](https://www.datanovia.com/en/lessons/anova-in-r/)
    -   [Repeated Measures ANOVA in
        R](https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/)
    -   [Mixed ANOVA in
        R](https://www.datanovia.com/en/lessons/mixed-anova-in-r/)
    -   [ANCOVA in R](https://www.datanovia.com/en/lessons/ancova-in-r/)
    -   [One-Way MANOVA in
        R](https://www.datanovia.com/en/lessons/one-way-manova-in-r/)
    -   [Kruskal-Wallis Test in
        R](https://www.datanovia.com/en/lessons/kruskal-wallis-test-in-r/)
    -   [Friedman Test in
        R](https://www.datanovia.com/en/lessons/friedman-test-in-r/)
