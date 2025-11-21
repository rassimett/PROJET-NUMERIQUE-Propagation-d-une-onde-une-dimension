function visualize_wave(filename)
    % --- Robust loader --------------------------------------------------
    printf("Chargement du fichier %s ...\n", filename);
    data = dlmread(filename, ",", 1, 0);     % skip header

    t = data(:,1);
    x = data(:,2);
    y = data(:,3);

    unique_t = unique(t);
    unique_x = unique(x);

    Nx = length(unique_x);
    Nt = length(unique_t);

    printf("Nx = %d, Nt = %d\n", Nx, Nt);

    % Reshape
    Y = reshape(y, Nx, Nt);

    % --- Figure ----------------------------------------------------------
    fig = figure("name", ["Wave Viewer : " filename], "numbertitle", "off");
    set(fig, "color", "w");

    % fix scaling
    ymin = min(Y(:));
    ymax = max(Y(:));

    % If too many frames, skip for speed
    max_frames = 500;     % adjust if needed
    step = max(1, floor(Nt / max_frames));

    printf("Frame step = %d (acceleration ON)\n", step);

    % --- Animation -------------------------------------------------------
    for k = 1:step:Nt

        if ~ishandle(fig)
            printf("Fenetre fermee -> animation stoppee.\n");
            return;
        end

        plot(unique_x, Y(:,k), "linewidth", 2);
        ylim([ymin ymax]);
        xlim([0 unique_x(end)]);
        grid on;

        title(sprintf("Propagation — t = %.4f s  (%d/%d)", ...
              unique_t(k), k, Nt));

        xlabel("Position x (m)");
        ylabel("Deplacement y (m)");

        drawnow;

    endfor

    printf("Animation terminee.\n");
endfunction
