function [pref, pval] = get_unit_preference(fnd, cond, opt)

def.epoch = 1;
def.t_range = [100 200];
def.require_significance = false;

opt = safeStructAssign(def, opt);

ucond = unique(cond(~isnan(cond(:))));
if length(ucond) ~= 2
    error('there should be 2 conditions to compare.');
end


FR = fnd.FR({opt.epoch, opt.t_range});

if ~isequal(size(FR), size(cond))
    error('The number of units, trials in cond does not match fnd.');
end

ncell = size(FR,1);

pref = nan(ncell, 1);
pval = nan(ncell, 1);

for n = 1:ncell
    FR1 = FR(n, cond(n,:) == ucond(1));
    FR2 = FR(n, cond(n,:) == ucond(2));
    FR1 = FR1(~isnan(FR1));
    FR2 = FR2(~isnan(FR2));
    if isempty(FR1) || isempty(FR2)
        warning('Only one condition exists in unit %d. No selectivity defined.', n);
    end
    if mean(FR1) > mean(FR2)
        pref(n) = ucond(1);
    elseif mean(FR1) < mean(FR2)
        pref(n) = ucond(2);
    end
    [~, pval(n)] = ttest2(FR1, FR2);
end

if opt.require_significance
    pref(pval > .05) = NaN;
end



