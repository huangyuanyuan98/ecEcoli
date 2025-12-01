function savePDF(fig, filename)
% savePDF Save a figure as PDF without resizing
%
% Usage:
%   savePDF(fig, 'crabtree.pdf')
%
% Inputs:
%   fig      - Figure handle, e.g., gcf
%   filename - Full path for saving the PDF

    if nargin < 1 || isempty(fig), fig = gcf; end
    if nargin < 2, error('A PDF filename must be provided.'); end

    fig.Units = 'inches';
    sz = fig.Position(3:4);

    fig.PaperUnits = 'inches';
    fig.PaperSize = sz;
    fig.PaperPosition = [0 0 sz];

    print(fig, filename, '-dpdf', '-painters');
end
