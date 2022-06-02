//
// Copyright (C) 2022 Armin Sobhani <arminms@gmail.com>
//
// The MIT License
//

#include <iostream>
#include <cstdio>
#include <iomanip>
#include <unordered_map>
#include <numeric>

#include <zlib.h>
#include <kseq.h>
#include <cxxopts.hpp>

#include <version.hpp>

KSEQ_INIT(gzFile, gzread)

int main(int argc, char* argv[])
{
    // turning sync off before any standard i/o occurs
    std::cout.sync_with_stdio(false);

    try
    {
        cxxopts::Options options("ngx", " (amstools) -- print contig statistics\n");
        options.custom_help(
            "[OPTION]... [FILE]...\n"
            "  ngx [OPTION]... --files-from=F\n\n"
            "Print the contiguity statistics (e.g. N50, L50) for each FILE.\n"
            "Both FastA and FastQ (optionally gzipped) files are supported.\n"
            "Print NG/LG variants if expected genome size is provided.\n\n"
            "With no FILE, or when FILE is -, read standard input.\n\n"
            "The options below may be used to select which statistics are "
            "printed,\nalways in the following order: #Seq, #Res, Min, Max,"
            " N(G)x..., L(G)x..., File.\n"
        );
        options.add_options()
        (   "g,genome-size"
        ,   "expected genome size\n"
            "  if G is provided then NGx/LGx values\n"
            "  will be computed"
        ,   cxxopts::value<size_t>()
        ,   "G"
        )
        (   "f,files-from"
        ,   "read input from the files specified by\n"
            "  names separated by newlines in file F\n"
            "  If F is - then read names from standard input"
        ,   cxxopts::value<std::string>()
        ,   "F"
        )
        (   "l,lx-values"
        ,   "print Lx along with Nx values"
        )
        (   "m,min"
        ,   "minimum contig length to be considered\n"
            "  every contig sequence of length shorter\n"
            "  than M will be discarded\n "
        ,   cxxopts::value<size_t>()
        ->  default_value("1")
        ,   "M"
        )
        (   "n,nx-values"
        ,   "Nx values to be printed (e.g. -n50,90 for N50\n"
            "  and N90)"
        ,   cxxopts::value<std::vector<size_t>>()
        ->  default_value("50")
        ,   "x..."
        )
        (   "s,sequence-lengths"
        ,   "print sequence lengths statistics"
        )
        (   "help"
        ,   "display this help and exit"
        )
        (   "version"
        ,   "output version information and exit"
        )
        (   "files"
        ,   "files"
        ,   cxxopts::value<std::vector<std::string>>()
        )
        ;

        options.parse_positional({"files"});
        options.positional_help("");

        auto result = options.parse(argc, argv);

        if (result.count("help"))
        {
            std::cout << options.program()
                      << options.help()
                      << std::endl;
            return 0;
        }

        if (result.count("version"))
        {
            std::cout << options.program()
                      << AMSTOOLS_VERSION
                      << std::endl;
            return 0;
        }

        if (result.count("files") && result.count("files-from"))
        {
            std::cerr << options.program() << ": "
                        << "file operands cannot be combined with --files-from"
                        << std::endl;
            return 1;
        }

        // making file list
        std::vector<std::string> file_stdin{ "-" };
        auto& files_in = result.count("files")
        ?   result["files"].as<std::vector<std::string>>()
        :   file_stdin;

        std::vector<std::string> files_from;
        if (result.count("files-from"))
        {
            auto& file = result["files-from"].as<std::string>();
            gzFile fp = file == "-"
            ?   gzdopen(fileno(stdin), "r")
            :   gzopen(file.c_str(), "r");
            if (nullptr == fp)
                std::cerr << options.program() << ": "
                            << "error reading "
                            << file
                            << std::endl;
            kstream_t* ks = ks_init(fp);
            kstring_t str = {0,0,0};
            while (ks_getuntil(ks, '\n', &str, 0) >= 0)
                files_from.emplace_back(str.s);
            ks_destroy(ks);
            gzclose(fp);
            free(str.s);
        }

        auto& files = result.count("files-from") ? files_from : files_in;
        auto& threshold = result["nx-values"].as<std::vector<size_t>>();

        // printing header
        std::cout << std::setw(11) << std::left << "#Seq"
                  << std::setw(11) << std::left << "#Res";
        if (result.count("sequence-lengths"))
        {
            std::cout << std::setw(11) << std::left << "Min"
                      << std::setw(11) << std::left << "Max";
        }
        for (size_t i = 0; i < threshold.size(); ++i)
        {
            std::string ng = result.count("genome-size") ? "NG" : "N";
            ng += std::to_string(threshold[i]);
            std::cout << std::setw(11) << std::left << ng;
        }
        if (result.count("lx-values"))
        {
            for (size_t i = 0; i < threshold.size(); ++i)
            {
                std::string lg = result.count("genome-size") ? "LG" : "L";
                lg += std::to_string(threshold[i]);
                std::cout << std::setw(11) << std::left << lg;
            }
        }
        std::cout << "File" << std::endl;

        std::vector<size_t> contig_length;
        for (const auto& file : files)
        {
            gzFile fp = file == "-"
            ?   gzdopen(fileno(stdin), "r")
            :   gzopen(file.c_str(), "r");
            if (nullptr == fp)
            {
                std::cerr << options.program() << ": "
                            << "error reading:\t"
                            << file
                            << std::endl;
                continue;
            }
            kseq_t* seq = kseq_init(fp);
            std::unordered_map<char, size_t> bp_counter(7);
            while (kseq_read(seq) >= 0)
                contig_length.push_back(seq->seq.l);
            kseq_destroy(seq);
            gzclose(fp);

            // 1. ordering contigs by their lengths from the longest to the
            //    shortest
            std::sort(
                contig_length.begin()
            ,   contig_length.end()
            ,   std::greater<size_t>()
            );

            // 2. set the number of contigs and residues
            auto min_length = result["min"].as<size_t>();
            auto n_contigs = contig_length.size();
            if (min_length > 1)
                for (size_t i = 0; i < contig_length.size(); ++i)
                    if (contig_length[i] < min_length)
                    {
                        n_contigs = i;
                        break;
                    }
            auto n_res = std::accumulate(
                contig_length.begin()
            ,   contig_length.begin() + n_contigs
            ,   0 );

            // 3. calculate the cutoff value by summing all contigs and
            //    multiplying by the threshold percentage
            std::vector<size_t> cutoff;
            cutoff.reserve(threshold.size());
            if (result.count("genome-size"))
            {
                auto genome_size = result["genome-size"].as<size_t>();
                for (auto x : threshold)
                    cutoff.push_back(genome_size * x / 100);
            }
            else
            {
                for (auto x : threshold)
                    cutoff.push_back(n_res * x / 100);
            }

            // 4. compute LN(G)x values
            std::vector<size_t> lgx_value(cutoff.size());
            std::vector<size_t> ngx_value(cutoff.size());
            for (size_t i = 0; i < cutoff.size(); ++i)
            {
                size_t sum{};
                for (size_t j = 0; j < n_contigs; ++j)
                {
                    lgx_value[i] = j + 1;
                    ngx_value[i] = contig_length[j];
                    sum += contig_length[j];
                    if (sum >= cutoff[i])
                        break;
                }
            }

            // printing values
            std::cout << std::setw(11) << std::left << n_contigs
                      << std::setw(11) << std::left << n_res;
            if (result.count("sequence-lengths"))
                std::cout << std::setw(11) << std::left
                          << contig_length[n_contigs - 1] // min
                          << std::setw(11) << std::left
                          << contig_length[0];  // max
            for (size_t i = 0; i < threshold.size(); ++i)
                std::cout << std::setw(11) << std::left
                          << ngx_value[i];
            if (result.count("lx-values"))
                for (size_t i = 0; i < threshold.size(); ++i)
                    std::cout << std::setw(11) << std::left
                            << lgx_value[i];
            std::cout << file << std::endl;
        }
    }
    catch(const cxxopts::OptionException& e)
    {
        std::cerr << "ngx: " << e.what() << std::endl;
        return 1;
    }
    catch (std::exception& e)
    {
        std::cerr << "ngx: " << e.what() << std::endl;
    }
}
