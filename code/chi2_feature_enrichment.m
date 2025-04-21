function [p_values, bh_corrected_p, chi2_stats, is_significant] = chi2_feature_enrichment(feature_counts, C, D, fdr_threshold)
    if nargin < 4
        fdr_threshold = 0.05; % Default FDR threshold
    end

    % Number of features
    m = size(feature_counts, 1);

    % Initialize output
    p_values = zeros(m, 1);
    chi2_stats = zeros(m, 1);

    % Loop over features and compute p-values
    for i = 1:m
        % Counts
        a = feature_counts(i, 1);  % present in A
        b = feature_counts(i, 2);  % present in B
        c = C - a;                 % absent in A
        d = D - b;                 % absent in B

        % Construct contingency table
        table = [a, b; c, d];

        % Perform Chi-squared test
        [~, p, stats] = chi2gof_from_table(table);

        % Store results
        p_values(i) = p;
        chi2_stats(i) = stats.chi2stat;
    end

    % Benjamini-Hochberg correction
    [bh_corrected_p, is_significant] = benjamini_hochberg(p_values, fdr_threshold);
end

function [h, p, stats] = chi2gof_from_table(contingency)
    % Compute expected frequencies
    row_totals = sum(contingency, 2);
    col_totals = sum(contingency, 1);
    total = sum(row_totals);

    expected = (row_totals * col_totals) / total;

    % Chi-squared statistic
    chi2stat = sum((contingency - expected).^2 ./ expected, 'all');

    % Degrees of freedom for 2x2 table = 1
    df = 1;

    % p-value
    p = 1 - chi2cdf(chi2stat, df);

    % Output
    h = p < 0.05;  % default alpha
    stats.chi2stat = chi2stat;
    stats.df = df;
end

function [bh_corrected_p, is_significant] = benjamini_hochberg(p_values, fdr)
    m = length(p_values);
    [sorted_p, sort_idx] = sort(p_values);
    
    % Adjusted p-values
    adjusted_p = zeros(m, 1);
    for i = m:-1:1
        if i == m
            adjusted_p(i) = sorted_p(i);
        else
            adjusted_p(i) = min(sorted_p(i) * m / i, adjusted_p(i+1));
        end
    end

    % Reorder to original order
    bh_corrected_p = zeros(m,1);
    bh_corrected_p(sort_idx) = adjusted_p;

    % Significance flag
    is_significant = bh_corrected_p < fdr;
end
