# A few headers displays
function print_header(Title)

    @printf(" \n\n %-15s\n", Title)
    @printf(
        "\n%-16s  %5s  %9s  %7s %7s  %5s  %5s  %5s  %6s  %s   %s\n",
        "Name",
        "nvar",
        "f",
        "‖∇f‖∞",
        "‖∇f‖₂",
        "#obj",
        "#grad",
        " H  ",
        "#hprod",
        "status",
        "time"
    )
end


function display_header_problems()
    @printf(
        "%-15s  %5s  %9s  %7s %7s  %5s  %5s  %6s  %s\n",
        "Name",
        "nvar",
        "f",
        "‖∇f‖∞",
        "‖∇f‖₂",
        "#obj",
        "#grad",
        "#hprod",
        "status"
    )
end

function display_header_iterations()
    @printf(" Iter    f          ||g||       λ        α       status \n",)
end

# A few displays for verbose iterations
function display_failure(iter, fnext, λ, α)
    @printf("%4d  %10.3e             %7.1e %7.1e    unsuccessful\n", iter, fnext, λ, α)
end

function display_v_success(iter, f, norm_g, λ, α)
    @printf("%4d  %10.3e %9.2e   %7.1e %7.1e Very successful\n", iter, f, norm_g, λ, α)
end

function display_success(iter, f, norm_g, λ, α)
    @printf("%4d  %10.3e %9.2e   %7.1e %7.1e      successful\n", iter, f, norm_g, λ, α)
end

function print_stats(prob, dim, f, gNorm, gnorm2, calls, status, timt)
    @printf(
        "%-16s  %5d  %9.2e  %7.1e  %7.1e  %5d  %5d %5d %6d  %s  %8.3f\n",
        prob,
        dim,
        f,
        gNorm,
        gnorm2,
        calls[1],
        calls[2],
        calls[3],
        calls[4],
        status,
        timt
    )
end
