//
// Copyright (C) 2022 Armin Sobhani <arminms@gmail.com>
//
// The MIT License
//

#include <iostream>
#include <cstdio>

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
        cxxopts::Options options("sc", AMSTOOLS_TITLE);
        options.custom_help(
            "[OPTION]... [FILE]...\n"
            "  sc [OPTION]... --files-from=F\n"
            "Print seqs and bps counts for each FILE, and total values if more than one\n"
            "FILE is specified. Both FastA and FastQ (optionally gzipped) files are supported.\n\n"
            "With no FILE, or when FILE is -, read standard input.\n\n"
            "The options below may be used to select which counts are printed, always in\n"
            "the following order: seqs, bps, maximum sequence length."
        );
        options.add_options()
        ("b,bps", "print the base pair counts")
        ("s,seqs", "print the sequence counts")
        ("m,max-seq-length", "print the longest sequence counts")
        ("f,files-from"
        ,   "read input from the files specified by\n"
            "  names separated by newlines in file F;\n"
            "  If F is - then read names from standard input"
        ,   cxxopts::value<std::string>(), "F" )
        ("help", "display this help and exit")
        ("version", "output version information and exit")
        ("files", "files", cxxopts::value<std::vector<std::string>>())
        ;

        options.parse_positional({"files"});
        options.positional_help("");

        auto result = options.parse(argc, argv);

        if (result.count("help"))
        {
            std::cout << options.help() << std::endl;
            return 0 ;
        }

        if (result.count("version"))
        {
            std::cout << AMSTOOLS_VERSION
                      << std::endl;
            return 0 ;
        }

        if (result.count("files") || result.count("files-from"))
        {
            if (result.count("files") && result.count("files-from"))
            {
                std::cerr << "sc: file operands cannot be combinded with --files-from"
                          << std::endl;
                return 1;
            }

            std::vector<std::string> files_from;
            if (result.count("files-from"))
            {
                auto& file = result["files-from"].as<std::string>();
                gzFile fp = file == "-"
                ?   gzdopen(fileno(stdin), "r")
                :   gzopen(file.c_str(), "r");
                if (nullptr == fp)
                    std::cerr << "sc: error reading " << file << std::endl;
                kstream_t* ks = ks_init(fp);
                kstring_t str = {0,0,0};
                while (ks_getuntil(ks, '\n', &str, 0) >= 0)
                    files_from.emplace_back(str.s);
                ks_destroy(ks);
                gzclose(fp);
                free(str.s);
            }

            auto& files = result.count("files")
            ?   result["files"].as<std::vector<std::string>>()
            :   files_from;
            size_t seqsn_total{}, bpsn_total{}, seqmax_total{};
            for (const auto& file : files)
            {
                size_t seqsn{}, bpsn{}, seqmax{};
                gzFile fp = file == "-"
                ?   gzdopen(fileno(stdin), "r")
                :   gzopen(file.c_str(), "r");
                if (nullptr == fp)
                {
                    std::cerr << "error reading:\t\t" << file << std::endl;
                    continue;
                }
                kseq_t* seq = kseq_init(fp);
                while (kseq_read(seq) >= 0)
                {
                    ++seqsn;
                    bpsn += seq->seq.l;
                    if (seq->seq.l > seqmax)
                        seqmax = seq->seq.l;
                }
                seqsn_total += seqsn;
                bpsn_total += bpsn;
                if (seqmax > seqmax_total)
                    seqmax_total = seqmax;
                kseq_destroy(seq);
                gzclose(fp);
                if (0 == result.count("seqs")
                &&  0 == result.count("bps")
                &&  0 == result.count("max-seq-length") )
                    std::cout << seqsn << '\t'
                              << bpsn  << '\t'
                              << file
                              << std::endl;
                else
                {
                    if (result.count("seqs"))
                        std::cout << seqsn << '\t';
                    if (result.count("bps"))
                        std::cout << bpsn << '\t';
                    if (result.count("max-seq-length"))
                        std::cout << seqmax << '\t';
                    std::cout << file << std::endl;
                }
            }
            if (files.size() > 1)
            {
                if (0 == result.count("seqs")
                &&  0 == result.count("bps")
                &&  0 == result.count("max-seq-length") )
                    std::cout << seqsn_total << '\t'
                              << bpsn_total
                              << "\ttotal"
                              << std::endl;
                else
                {
                    if (result.count("seqs"))
                        std::cout << seqsn_total << '\t';
                    if (result.count("bps"))
                        std::cout << bpsn_total << '\t';
                    if (result.count("max-seq-length"))
                        std::cout << seqmax_total << '\t';
                    std::cout << "total" << std::endl;
                }
            }
        }
        else
        {
            size_t seqsn{}, bpsn{}, seqmax{};
            gzFile fp = gzdopen(fileno(stdin), "r");
            kseq_t *seq = kseq_init(fp);
            while (kseq_read(seq) >= 0)
            {
                ++seqsn;
                bpsn += seq->seq.l;
                if (seq->seq.l > seqmax)
                    seqmax = seq->seq.l;
            }
            kseq_destroy(seq);
            gzclose(fp);
            if (0 == result.count("seqs")
            &&  0 == result.count("bps")
            &&  0 == result.count("max-seq-length") )
                std::cout << seqsn << '\t' << bpsn << std::endl;
            else
            {
                if (result.count("seqs"))
                    std::cout << seqsn << '\t';
                if (result.count("bps"))
                    std::cout << bpsn << '\t';
                if (result.count("max-seq-length"))
                    std::cout << seqmax << '\t';
                std::cout << std::endl;
            }
        }
    }
    catch(const cxxopts::OptionException& e)
    {
        std::cerr << "sc: " << e.what() << std::endl;
        return 1;
    }
}