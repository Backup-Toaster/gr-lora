/* -*- c++ -*- */
/*
 * Copyright 2018 Pieter Robyns.
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include <gnuradio/expj.h>
#include <lora/utilities.h>
#include "encoder_impl.h"
#include "tables.h"
#include <bitset> //LT Edit.

namespace gr {
  namespace lora {

    encoder::sptr
    encoder::make() {
      return gnuradio::get_initial_sptr
        (new encoder_impl());
    }

    /*
     * The private constructor
     */
    encoder_impl::encoder_impl() : gr::sync_block("encoder",
              gr::io_signature::make(0, 0, 0),
              gr::io_signature::make(1, 1, sizeof(gr_complex))) {
        // Bind ports
        message_port_register_in(pmt::mp("in"));
        set_msg_handler(pmt::mp("in"), boost::bind(&encoder_impl::handle_loratap, this, _1));

        set_output_multiple(pow(2,16));

        // Initialize variables
        d_osr = 8;
        d_samples_per_second = 125000*d_osr;
        d_num_preamble_symbols = 8;
        d_bw = 125000;
        d_explicit = true;
        d_reduced_rate = false;
        d_chirp_phi0 = -M_PI;

        d_dt = 1.0f / d_samples_per_second;
        d_sample_buffer.reserve(d_samples_per_second); // Allocate for one second of samples

        // Setup hamming coding schemes
        fec_scheme fs = LIQUID_FEC_HAMMING84;
        d_h48_fec = fec_create(fs, NULL);

        // Setup chirp lookup tables TODO unused now
        for(uint8_t sf = 6; sf <= 12; sf++) {
            const uint32_t chips_per_symbol = pow(2, sf);
            const double symbols_per_second = (double)d_bw / chips_per_symbol;
            const double samples_per_symbol = d_samples_per_second / symbols_per_second;
            const double T = 0.5 * d_bw * symbols_per_second;
            const double f0 = -(d_bw / 2.0);
            const double pre_dir = 2.0 * M_PI;
            double t;

            std::vector<gr_complex> chirp(samples_per_symbol*2);
            for (uint32_t i = 0u; i < samples_per_symbol; i++) {
                t = d_dt * i;
                gr_complex sample = gr_expj(pre_dir * t * (f0 + T * t));
                chirp[i] = sample;
                chirp[i+samples_per_symbol] = sample;
            }

            std::cout << "SF " << (int)sf << " has " << samples_per_symbol << " samples" << std::endl;

            d_chirps[sf] = chirp; // Copy vector metadata to chirps hashmap. Note that vector internally allocates on heap.
        }
    }

    void encoder_impl::transmit_chirp(bool up, uint8_t sf, uint16_t symbol, bool quarter = false) {
        const uint32_t chips_per_symbol = pow(2, sf);
        const double symbols_per_second = (double)d_bw / chips_per_symbol;
        const uint32_t samples_per_symbol = d_samples_per_second / symbols_per_second;
        const double T = 0.5 * d_bw * symbols_per_second;
        const double f0 = -(d_bw / 2.0);
        double pre_dir = 2.0 * M_PI;
        double t = 0;
        // std::vector<gr_complex> chirp_0(2*samples_per_symbol);
        // chirp_0 = d_chirps[sf];
        // std::cout << "chirp_0.size(): " << chirp_0.size() << std::endl;

        if(!up)
            pre_dir *= -1;

        std::vector<gr_complex> chirp(samples_per_symbol);
        double phase = 0;
        // std::cout << "d_osr*symbol + 1u: " << d_osr*symbol + 1u << std::endl;
        // std::cout << "(2*samples_per_symbol - 1) - d_osr*symbol - 1u: " << (2*samples_per_symbol - 1) - d_osr*symbol - 1u << std::endl;
        // std::cout << "2*samples_per_symbol - d_osr*symbol: " << 2*samples_per_symbol - d_osr*symbol << std::endl;
        // if (!up) {
        //     for (uint32_t i = 0u; i < samples_per_symbol; i++) {
        //         chirp[i] = chirp_0[d_osr*symbol + i];
        //         // chirp[i] = chirp_0[(2*samples_per_symbol - 1) - d_osr*symbol - i];
        //     }
        // } 
        // else {
        //     for (uint32_t i = 0u; i < samples_per_symbol; i ++) {
        //         chirp[i] = chirp_0[(2u*samples_per_symbol - 1u) - d_osr*symbol - i];
        //     }
        // }
        

        uint32_t offset = 0u;
        for (uint32_t i = 0u; i < samples_per_symbol; i++) { // Original
        // for (uint32_t i = 1u; i <= samples_per_symbol; i++) { // LT Edit - Shifted i up by one.
        // t = d_dt * ((i + (d_osr * symbol)) % samples_per_symbol); // Original
        t = d_dt * ((i + (d_osr * (symbol))) % samples_per_symbol); //LT Edit
        phase = d_chirp_phi0 + (pre_dir * t * (f0 + T * t));
        chirp[i] = gr_expj(phase); //Original
        chirp[i - offset] = gr_expj(phase); // LT Edit
        // std::cout << "i: " << i << ", t: " << t << ", phase: " << phase << ", chirp[i]: " << chirp[i - offset] << std::endl;
        }

        // Add chirp to buffer
        if(quarter)
            d_sample_buffer.insert(d_sample_buffer.end(), chirp.begin(), chirp.begin() + chirp.size() / 4);
        else
            d_sample_buffer.insert(d_sample_buffer.end(), chirp.begin(), chirp.end());

        // Set phase
        d_chirp_phi0 = fmod(phase,(2 * M_PI));
    }

    /*
     * Our virtual destructor.
     */
    encoder_impl::~encoder_impl() {
        fec_destroy(d_h48_fec);
    }

    void encoder_impl::handle_loratap(pmt::pmt_t msg) {
        gr::thread::scoped_lock guard(d_mutex); // Auto lock & unlock (RAII)

        d_packets.push_back(msg);
    }

    bool encoder_impl::parse_packet_conf(loraconf_t& conf, uint8_t* packet, uint32_t packet_len) {
        if(packet_len <= sizeof(loraconf_t)) {
            return false;
        }

        memcpy(&conf.tap, packet, sizeof(loratap_header_t));
        memcpy(&conf.phy, packet + sizeof(loratap_header_t), sizeof(loraphy_header_t));

        return true;
    }

    /**
     * Get packets from queue and initialize transmission.
     */
    int encoder_impl::work(int noutput_items, gr_vector_const_void_star &input_items, gr_vector_void_star &output_items) {
        gr_complex* out = (gr_complex*)output_items[0];

        // Temporary         ve  pa      le              fq  bw  sf  pr  mr  cr  sn  sy  H1  H1  H1
        // char test_pkt[] = "\x00\x00\x12\x00\x00\xa1\xbc\x33\x01\x08\x00\x00\x00\x00\x12\x17\x94\xa0\x00\x01\x02\x03\x04\x05\x06\x07\x08\x09\x10\x11\x12\x13\x14\x15\x16\x17\x18\x19\x20\x21\x22\xb8\x73";
        // char test_pkt[] = "\x00\x00\x12\x00\x00\xa1\xbc\x33\x01\x07\x00\x00\x00\x00\x00\x01\x80\xa0\x0f";
        // char test_pkt[] = "\x00\x00\x12\x00\x00\xa1\xbc\x33\x01\x0a\x00\x00\x00\x00\x00\x04\x80\xa0\xde\xad\xbe\xef";
        char test_pkt[] = "\x00\x00\x12\x00\x00\xa1\xbc\x33\x01\x0a\x00\x00\x00\x00\x00\x12\x91\x60\x00\x00\x00\x00\x01\x50\x45\x52\x08\x09\x0a\x0b\x0c\x0d\x0e\x0f\x10\x11\xd4\xcc";
        // char test_pkt[] = "\x00\x00\x12\x00\x00\xa1\xbc\x33\x01\x08\x00\x00\x00\x00\x00\x08\x80\xa0\x0f\x0f\x0f\x0f\x0f\x0f\x0f\x0f";
        // char test_pkt[] = "\x00\x00\x12\x00\x00\xa1\xbc\x33\x01\x07\x00\x00\x00\x00\x12\x02\x91\xa0\x00\xff";
        // char test_pkt[] = "\x00\x00\x12\x00\x00\xa1\xbc\x33\x01\x07\x00\x00\x00\x00\x12\x17\x91\xa0\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00";
        // char test_pkt[] = "\x00\x00\x12\x00\x00\xa1\xbc\x33\x01\x07\x00\x00\x00\x00\x12\x17\x91\xa0\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02\x02";

        loraconf_t conf;
        memset(&conf, 0, sizeof(loraconf_t));

        if(!parse_packet_conf(conf, (uint8_t*)test_pkt, sizeof(test_pkt))) {
            std::cerr << "Malformed LoRa packet received" << std::endl;
            exit(1);
        }

        print_loraconf(conf);
        transmit_packet(conf, (uint8_t*)(test_pkt + sizeof(loratap_header_t)));

        if(d_sample_buffer.size() >= noutput_items) {
            // Get noutput_items from buffer
            memcpy(out, &d_sample_buffer[0], noutput_items * sizeof(gr_complex));
            d_sample_buffer.erase(d_sample_buffer.begin(), d_sample_buffer.begin() + noutput_items);

            return noutput_items;
        } else {
            return 0; // Wait for more symbols
        }
        // -------------------------------------------------------------------------------------

        // Get packet from queue
        gr::thread::scoped_lock guard(d_mutex);
        if(d_packets.size() > 0) {
            pmt::pmt_t packet = d_packets.front();

            // Process
            std::cout << "Processing packet" << std::endl;
            //uint8_t* packet = (uint8_t*)pmt::blob_data(msg);
            //size_t packet_length = pmt::blob_length(msg);

            // Delete from queue
            d_packets.erase(d_packets.begin());
        }

        // Tell runtime system how many output items we produced.
        return 0;
    }

    void encoder_impl::shuffle(uint8_t *data, uint32_t data_len, const uint8_t *shuffle_pattern) {
        for (uint32_t i = 0u; i < data_len; i++) {
            uint8_t result = 0u;

            for (uint32_t j = 0u; j < 8; j++) {
                result |= !!(data[i] & (1u << shuffle_pattern[j])) << j;
            }

            data[i] = result;
        }
    }

    uint32_t encoder_impl::interleave_block(uint16_t *symbols, uint8_t* data, uint8_t sf, uint8_t cr, bool reduced_rate) {
        if(reduced_rate)
            sf -= 2;

        // Determine symbols for this block
        for(uint8_t symbol = 0; symbol < 4+cr; symbol++) {
            for(int8_t bit = sf-1; bit >= 0; bit--) {
                symbols[symbol] |= ((data[bit] >> symbol) & 0x01) << bit;
            // std::cout << "Loop Output  " << std::bitset<16> (symbols[symbol]) << std::endl;
            }

            int32_t to_rotate = sf-symbol;
            if(to_rotate < 0)
                to_rotate += sf;

            symbols[symbol] = gr::lora::rotl(symbols[symbol], to_rotate, sf);
        }

        // Rotate to interleave
        std::vector<uint16_t> symbols_v(symbols, symbols + (4+cr));
        print_interleave_matrix(std::cout, symbols_v, sf);
        print_vector(std::cout, symbols, "Chips", 4+cr, sf);

        // Determine bins
        std::cout << "Bins (Original; Offset): ";

        // Original
        for(uint8_t symbol = 0; symbol < 4+cr; symbol++) {
            symbols[symbol] = gray_decode(symbols[symbol]);
            if(reduced_rate)
                symbols[symbol] <<= 2;
            symbols[symbol] += 1; //LT Addition
            std::cout << "(" << (int)symbols[symbol] - 1 << "; " << (int)symbols[symbol] << "), ";
        }
        std::cout << std::endl;

        return sf;
    }

    void encoder_impl::nibble_swap(uint8_t* encoded, uint32_t length) {
        for(uint32_t i = 0; i+1 < length; i += 2) {
            uint8_t tmp = encoded[i];
            encoded[i] = encoded[i+1];
            encoded[i+1] = tmp;
        }
    }

    void encoder_impl::transmit_packet(loraconf_t& conf, uint8_t* packet) { // TODO: clean up
        uint32_t header_length = sizeof(loraphy_header_t);
        uint32_t payload_length = conf.phy.length + MAC_CRC_SIZE * conf.phy.has_mac_crc;
        uint32_t packet_length = header_length + payload_length;

        uint32_t bytes_needed_header_words = header_length * 2; // 2 = (4 + 4)/4, which is true since header coding always is (8,4).
        uint32_t bytes_needed_payload_words = std::ceil(payload_length * (4 + conf.phy.cr) / 4) - (conf.tap.channel.sf - 7);
        // Since SF - 2 codewords get packed into the header symbols, but only 5 header codewords are generated and used, SF-2-5 payload codewords are not included in the payload transmission.
        std::cout << "bytes_needed_payload_words: " << bytes_needed_payload_words << std::endl;
        // Small algorithm to find the next multiple of SF greater than or equal to the number of required bytes for payload words.
        uint32_t N = bytes_needed_payload_words;
        uint32_t K = conf.tap.channel.sf;
        uint32_t bytes_needed_payload_words_padded =  (((N + K) % K) == 0) ? N : (N + K) - (N + K) % K;
        std::cout << "bytes_needed_payload_words_padded: " << bytes_needed_payload_words_padded << std::endl;

        uint32_t num_bytes = bytes_needed_header_words - 1 + bytes_needed_payload_words_padded + (conf.tap.channel.sf - 7); // Minus one as the last codeword of the header is overwritten.
        // We add back in that (SF - 7) value to make sure we are counting all our bytes.
        std::cout << "num_bytes: " << num_bytes << std::endl;

        uint32_t num_symbols_header = 8; //Equal to cr + 4, the output length of the interleaving function, when header cr is always 4.
        uint32_t num_symbols_payload = bytes_needed_payload_words_padded / conf.tap.channel.sf * (4 + conf.phy.cr); // Comes from fact that sf payload codewords are read in to produce 4 + cr symbols.

        uint32_t num_symbols = num_symbols_header + num_symbols_payload;
        std::cout << "num_symbols: " << num_symbols << std::endl;
        // uint32_t packet_length = conf.phy.length + sizeof(loraphy_header_t) + conf.phy.has_mac_crc * MAC_CRC_SIZE; // LT: Will look to see if CRC is included. 
        // uint32_t bits_needed_packet = 8 * packet_length;
        // uint32_t bits_needed_words = bits_needed_packet*(conf.phy.cr + 4)/conf.phy.cr; //When CR != 4 this will need to be split out for header and payload.
        // uint32_t num_symbols = std::ceil(bits_needed_words / conf.tap.channel.sf);
        // uint32_t num_bytes = bits_needed_words / 8;

        uint8_t encoded[num_bytes];
        memset(encoded, 0x00, num_bytes * sizeof(uint8_t)); // Set all symbols to zero so that zero symbols are transmitted if we don't use all the space.
        // uint32_t num_symbols = num_bytes * ((4.0+conf.phy.cr) / conf.tap.channel.sf) + 0.5;
        uint32_t encoded_offset = 0;
        uint32_t packet_offset = 0;

        // Add preamble symbols to queue
        for(uint32_t i = 0; i < d_num_preamble_symbols; i++) {
            transmit_chirp(true, conf.tap.channel.sf, 0);
        }

        // Add sync words to queue
        uint8_t sync_word = 0x00; // LT
        //uint8_t sync_word = 0x12; // Original
        uint32_t sync_offset_1 = ((sync_word & 0xf0) >> 4) * pow(2, conf.tap.channel.sf) * d_osr / 32;
        uint32_t sync_offset_2 = (sync_word & 0x0f) * pow(2, conf.tap.channel.sf) * d_osr / 32;
        transmit_chirp(true, conf.tap.channel.sf, sync_offset_1);
        transmit_chirp(true, conf.tap.channel.sf, sync_offset_2);

        // Add SFD to queue
        transmit_chirp(false, conf.tap.channel.sf, 0);
        transmit_chirp(false, conf.tap.channel.sf, 0);
        transmit_chirp(false, conf.tap.channel.sf, 0, true);

        // If explicit header, add one block (= SF codewords) to queue in reduced rate mode (and always 4/8)
        if(d_explicit) {
            fec_encode(d_h48_fec, 3, packet, encoded); // Header is always 4/8
            packet_offset = 3;
            // encoded_offset = 6; //LT
            encoded_offset = 5; // Original
        }

        // Add remaining blocks to queue
        print_vector(std::cout, packet, "Packet", packet_length, 8);

        fec_encode(d_h48_fec, packet_length - packet_offset, packet+packet_offset, encoded+encoded_offset); // TODO: change to appropriate scheme
        print_vector(std::cout, encoded, "Encoded", num_bytes, 8);
        nibble_swap(encoded+encoded_offset, num_bytes-encoded_offset);
        print_vector(std::cout, encoded, "Nibble Swap", num_bytes, 8);

        whiten(encoded+encoded_offset, gr::lora::prng_payload, num_bytes-encoded_offset);
        print_vector(std::cout, encoded, "Whitened", num_bytes, 8);

        const uint8_t shuffle_pattern[] = {1, 2, 3, 5, 4, 0, 6, 7};
        shuffle(encoded, num_bytes, shuffle_pattern);
        print_vector(std::cout, encoded, "Shuffled", num_bytes, 8);

        // Interleaving
        uint16_t symbols[num_symbols];
        memset(symbols, 0x00, num_symbols * sizeof(uint16_t));
        uint32_t symbols_done = 0;
        uint32_t interleave_offset = 0;

        if(d_explicit) {
            interleave_offset += interleave_block(symbols, encoded, conf.tap.channel.sf, 4, true);
            symbols_done += 8;
        }

        print_vector(std::cout, symbols, "Interleaved Header Only", symbols_done, conf.tap.channel.sf);

        while(interleave_offset+conf.tap.channel.sf <= num_bytes) { // TODO: needs to be exact number of bytes
            interleave_offset += interleave_block(symbols+symbols_done, encoded+interleave_offset, conf.tap.channel.sf, conf.phy.cr, d_reduced_rate);
            symbols_done += 4 + conf.phy.cr;
        }

        print_vector(std::cout, symbols, "Interleaved Symbols", symbols_done, conf.tap.channel.sf);
        
        // Transmission
        for(uint32_t i = 0; i < num_symbols; i++) {
            transmit_chirp(true, conf.tap.channel.sf, symbols[i]);
        }
    }

  } /* namespace lora */
} /* namespace gr */
