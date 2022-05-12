//
// Copyright (C) 2022 Armin Sobhani <arminms@gmail.com>
//
// The MIT License
//

#include <iostream>
#include <cstdio>
#include <iomanip>
#include <unordered_map>

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
        cxxopts::Options options("acgt", " (amstools) -- print residue statistics\n");
        options.custom_help(
            "[OPTION]... [FILE]...\n"
            "  acgt [OPTION]... --files-from=F\n\n"
            "Print residue statistics and optionally GC and AT contents for each"
            " FILE.\nBoth FastA and FastQ (optionally gzipped) files are "
            "supported.\n\nWith no FILE, or when FILE is -, read standard "
            "input.\n\nThe options below may be used to select which statistics "
            "are printed,\nalways in the following order: #seq, #res, residue "
            "statistics, AT-Content,\nGC-Content."
        );
        options.add_options()
        (   "a,AT-Content"
        ,   "print AT-Content percent"
        )
        (   "g,GC-Content"
        ,   "print GC-Content percent"
        )
        (   "f,files-from"
        ,   "read input from the files specified by\n"
            "  names separated by newlines in file F;\n"
            "  If F is - then read names from standard input"
        ,   cxxopts::value<std::string>()
        ,   "F"
        )
        (   "r,residues"
        ,   "list of characters to count as residues;\n"
            "  If R is 'all' then count all characters"
        ,   cxxopts::value<std::string>()
        ->  default_value("ACGT")
        ->  implicit_value("ACGT")
        ,   "R"
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
            return 0 ;
        }

        if (result.count("version"))
        {
            std::cout << options.program()
                      << AMSTOOLS_VERSION
                      << std::endl;
            return 0 ;
        }

        if (result.count("files") && result.count("files-from"))
        {
            std::cerr << options.program() << ": "
                        << "file operands cannot be combined with --files-from"
                        << std::endl;
            return 1;
        }

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

        for (const auto& file : files)
        {
            size_t seqsn{}, bpsn{};
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
            {
                ++seqsn;
                bpsn += seq->seq.l;
                for (size_t i = 0; i < seq->seq.l; ++i)
                    bp_counter[seq->seq.s[i]]++;
            }
            kseq_destroy(seq);
            gzclose(fp);

            // printing headers
            auto& residues = result["residues"].as<std::string>();
            std::cout << std::setw(10) << std::left << "#Seq" << ' '
                        << std::setw(10) << std::left << "#Res" << ' ';
            if ( result.count("residues")
            || (0 == result.count("GC-Content")
            &&  0 == result.count("AT-Content")
            &&  0 == result.count("residues") ) )
            {
                if (residues == "all")
                {
                    for (auto res : bp_counter)
                        std::cout << '#' << std::setw(9) << std::left
                                    << res.first << ' ';
                    for (auto res : bp_counter)
                        std::cout << '%' << std::setw(5) << std::left
                                    << res.first << ' ';
                }
                else
                {
                    for (size_t i = 0; i < residues.size(); ++i)
                        std::cout << '#' << std::setw(9) << std::left
                                    << residues[i] << ' ';
                    for (size_t i = 0; i < residues.size(); ++i)
                        std::cout << '%' << std::setw(5) << std::left
                                    << residues[i] << ' ';
                }
            }
            if (result.count("AT-Content"))
                std::cout << std::setw(6) << std::left << "%AT" << ' ';
            if (result.count("GC-Content"))
                std::cout << std::setw(6) << std::left << "%GC" << ' ';
            std::cout << "File" << std::endl;

            // printing values
            std::cout << std::setw(10) << std::left << seqsn << ' '
                        << std::setw(10) << std::left << bpsn << ' ';
            if ( result.count("residues")
            || (0 == result.count("GC-Content")
            &&  0 == result.count("AT-Content")
            &&  0 == result.count("residues") ) )
            {
                if (residues == "all")
                {
                    for (auto res : bp_counter)
                        std::cout << std::setw(10)
                                    << std::left
                                    << res.second
                                    << ' '
                                    ;
                    for (auto res : bp_counter)
                        std::cout << std::fixed
                                    << std::setw(5)
                                    << std::setprecision(2)
                                    << double(res.second)/bpsn*100
                                    << "% "
                                    ;
                }
                else
                {
                    for (size_t i = 0; i < residues.size(); ++i)
                        std::cout << std::setw(10) << std::left
                                    << bp_counter[residues[i]] << ' ';
                    for (size_t i = 0; i < residues.size(); ++i)
                        std::cout << std::fixed
                                    << std::setw(5)
                                    << std::setprecision(2)
                                    << double(bp_counter[residues[i]])/bpsn*100
                                    << "% "
                                    ;
                }
            }
            if (result.count("AT-Content"))
                std::cout  << std::fixed
                            << std::setw(5)
                            << std::setprecision(2)
                            <<     (double(bp_counter['A'])
                                +   double(bp_counter['T']))
                                /  (double(bp_counter['A'])
                                +   double(bp_counter['T'])
                                +   double(bp_counter['G'])
                                +   double(bp_counter['C']))
                                *   100
                            <<   "% "
                            ;
            if (result.count("GC-Content"))
                std::cout  << std::fixed
                            << std::setw(5)
                            << std::setprecision(2)
                            <<     (double(bp_counter['G'])
                                +   double(bp_counter['C']))
                                /  (double(bp_counter['A'])
                                +   double(bp_counter['T'])
                                +   double(bp_counter['G'])
                                +   double(bp_counter['C']))
                                *   100
                            <<   "% "
                            ;
            std::cout << file << std::endl;
        }
    }
    catch(const cxxopts::OptionException& e)
    {
        std::cerr << "acgt: " << e.what() << std::endl;
        return 1;
    }
}
