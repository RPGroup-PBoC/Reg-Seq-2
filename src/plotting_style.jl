using CairoMakie

function default_makie!()
    if ~isfile(assetpath("fonts", "Lucida-sans-Unicode-Regular.ttf"))
        #@warn "Lucida sans Unicode Regular font not added to Makie Fonts. Add to `~/.julia/packages/Makie/gQOQF/assets/fonts/`. Defaulting to NotoSans."
        Font = assetpath("fonts", "NotoSans-Regular.tff")
    else
        Font = assetpath("fonts", "Lucida-Sans-Unicode-Regular.ttf")
    end
    
    # Seaborn colorblind
    colors = ["#0173b2", "#de8f05", "#029e73", "#d55e00", "#cc78bc", "#ca9161", "#fbafe4", "#949494", "#ece133", "#56b4e9"]

    theme = Theme(
        Axis = (
            backgroundcolor = "#E3DCD0",
 
            # Font sizes
            titlesize=12,
            xlabelsize=12,
            ylabelsize=12,
            xticklabelsize=9,
            yticklabelsize=9,

            # Font styles
            titlefont=Font,
            xticklabelfont=Font,
            yticklabelfont=Font,
            xlabelfont=Font,
            ylabelfont=Font,

            # Grid
            xgridwidth=1.25,
            ygridwidth=1.25,
            xgridcolor="white",
            ygridcolor="white",
            xminorgridcolor="white",
            yminorgridcolor="white",
            xminorgridvisible=true,
            xminorgridwidth=1.,
            yminorgridvisible=true,
            yminorgridwidth=1,

            # Box
            rightspinevisible=false,
            topspinevisible=false,

            # Colorscheme
            palette = (color = colors,)

        ),
        Legend = (
            titlesize=8,
            labelsize=8,
            bgcolor="#E3DCD0",
            rowgap=-5,
            labelfont=Font

        ),
        backgroundcolor="white",
        linewidth=1.25,

    )
    set_theme!(theme)
end