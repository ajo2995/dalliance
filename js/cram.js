/* -*- mode: javascript; c-basic-offset: 4; indent-tabs-mode: nil -*- */

// 
// Dalliance Genome Explorer
// (c) Thomas Down 2006-2011
//
// cram.js: indexed binary alignments
//

"use strict";

if (typeof(require) !== 'undefined') {
    var spans = require('./spans');
    var Range = spans.Range;
    var union = spans.union;
    var intersection = spans.intersection;

    var bin = require('./bin');
    var readInt = bin.readInt;
    var readShort = bin.readShort;
    var readByte = bin.readByte;
    var readInt64 = bin.readInt64;
    var readFloat = bin.readFloat;

    var lh3utils = require('./lh3utils');
    var readVob = lh3utils.readVob;
    var unbgzf = lh3utils.unbgzf;
    var reg2bins = lh3utils.reg2bins;
    var Chunk = lh3utils.Chunk;

    var jszlib = require('jszlib');
}


var CRAM_MAGIC = 0x14d4142;
var CRAI_MAGIC = 0x88B1F;

var CramFlags = {
    MULTIPLE_SEGMENTS:       0x1,
    ALL_SEGMENTS_ALIGN:      0x2,
    SEGMENT_UNMAPPED:        0x4,
    NEXT_SEGMENT_UNMAPPED:   0x8,
    REVERSE_COMPLEMENT:      0x10,
    NEXT_REVERSE_COMPLEMENT: 0x20,
    FIRST_SEGMENT:           0x40,
    LAST_SEGMENT:            0x80,
    SECONDARY_ALIGNMENT:     0x100,
    QC_FAIL:                 0x200,
    DUPLICATE:               0x400,
    SUPPLEMENTARY:           0x800
};

function CramFile() {
}


// Calculate the length (in bytes) of the CRAI ref starting at offset.
// Returns {nbin, length, minBlockIndex}
function _getCraiRefLength(uncba, offset) {
    var p = offset;
    var nbin = readInt(uncba, p); p += 4;
    for (var b = 0; b < nbin; ++b) {
        var bin = readInt(uncba, p);
        var nchnk = readInt(uncba, p+4);
        p += 8 + (nchnk * 16);
    }
    var nintv = readInt(uncba, p); p += 4;

    var minBlockIndex = 1000000000;
    var q = p;
    for (var i = 0; i < nintv; ++i) {
        var v = readVob(uncba, q); q += 8;
        if (v) {
            var bi = v.block;
            if (v.offset > 0)
                bi += 65536;

            if (bi < minBlockIndex)
                minBlockIndex = bi;
            break;
        }
    }
    p += (nintv * 8);

    return {
        minBlockIndex: minBlockIndex,
        nbin: nbin,
        length: p - offset
    };
}


function makeCram(data, crai, indexChunks, callback, attempted) {
    // Do an initial probe on the CRAM file to catch any mixed-content errors.
    data.slice(0, 10).fetch(function(header) {
        if (header) {
            return makeCram2(data, crai, indexChunks, callback, attempted);
        } else {
            return callback(null, "Couldn't access CRAM.");
        }
    }, {timeout: 5000});
}

function makeCram2(data, crai, indexChunks, callback, attempted) {
    var cram = new CramFile();
    cram.data = data;
    cram.crai = crai;
    cram.indexChunks = indexChunks;

    var minBlockIndex = cram.indexChunks ? cram.indexChunks.minBlockIndex : 1000000000;

    // Fills out cram.chrToIndex and cram.indexToChr based on the first few bytes of the CRAM.
    function parseCramHeader(r) {
        if (!r) {
            return callback(null, "Couldn't access CRAM");
        }

        var unc = unbgzf(r, r.byteLength);
        var uncba = new Uint8Array(unc);

        var magic = readInt(uncba, 0);
        if (magic != CRAM_MAGIC) {
            return callback(null, "Not a CRAM file, magic=0x" + magic.toString(16));
        }
        var headLen = readInt(uncba, 4);
        var header = '';
        for (var i = 0; i < headLen; ++i) {
            header += String.fromCharCode(uncba[i + 8]);
        }

        var nRef = readInt(uncba, headLen + 8);
        var p = headLen + 12;

        cram.chrToIndex = {};
        cram.indexToChr = [];
        for (var i = 0; i < nRef; ++i) {
            var lName = readInt(uncba, p);
            var name = '';
            for (var j = 0; j < lName-1; ++j) {
                name += String.fromCharCode(uncba[p + 4 + j]);
            }
            var lRef = readInt(uncba, p + lName + 4);
            cram.chrToIndex[name] = i;
            if (name.indexOf('chr') == 0) {
                cram.chrToIndex[name.substring(3)] = i;
            } else {
                cram.chrToIndex['chr' + name] = i;
            }
            cram.indexToChr.push(name);

            p = p + 8 + lName;
        }

        if (cram.indices) {
            return callback(cram);
        }
    }

    function parseCrai(header) {
        if (!header) {
            return "Couldn't access CRAI";
        }

        var uncba = new Uint8Array(header);

        var craiMagic = readInt(uncba, 0);
        if (craiMagic != CRAI_MAGIC) {
            return callback(null, 'Not a CRAI file, magic=0x' + craiMagic.toString(16));
        }

        var offset = 10;
        var unc = jszlib.inflateBuffer(header,offset,header.byteLength - offset);
        var s = ''
        var resultBB = new Uint8Array(unc);
        for (var i = 0; i < resultBB.length; ++i) {
          s += String.fromCharCode(resultBB[i]);
        }
        console.log(s);

        cram.indices = [];

        var p = 8;
        for (var ref = 0; ref < nref; ++ref) {
            var blockStart = p;
            var o = _getCraiRefLength(uncba, blockStart);
            p += o.length;

            minBlockIndex = Math.min(o.minBlockIndex, minBlockIndex);

            var nbin = o.nbin;

            if (nbin > 0) {
                cram.indices[ref] = new Uint8Array(header, blockStart, p - blockStart);
            }
        }

        return true;
    }

    if (!cram.indexChunks) {
        cram.crai.fetch(function(header) {   // Do we really need to fetch the whole thing? :-(
            var result = parseCrai(header);
            if (result !== true) {
                if (cram.crai.url && typeof(attempted) === "undefined") {
                    // Already attempted x.cram.crai not there so now trying x.crai
                    cram.crai.url = cram.data.url.replace(new RegExp('.cram$'), '.crai');
                    
                     // True lets us know we are making a second attempt
                    makeCram2(data, cram.crai, indexChunks, callback, true);
                }
                else {
                    // We've attempted x.cram.crai & x.crai and nothing worked
                    callback(null, result);
                }
            } else {
              cram.data.slice(0, minBlockIndex).fetch(parseCramHeader);
            }
        });   // Timeout on first request to catch Chrome mixed-content error.
    } else {
        var chunks = cram.indexChunks.chunks;
        cram.indices = []
        for (var i = 0; i < chunks.length; i++) {
           cram.indices[i] = null;  // To be filled out lazily as needed
        }
        cram.data.slice(0, minBlockIndex).fetch(parseCramHeader);
    }
}



CramFile.prototype.blocksForRange = function(refId, min, max) {
    var index = this.indices[refId];
    if (!index) {
        return [];
    }

    var intBinsL = reg2bins(min, max);
    var intBins = [];
    for (var i = 0; i < intBinsL.length; ++i) {
        intBins[intBinsL[i]] = true;
    }
    var leafChunks = [], otherChunks = [];

    var nbin = readInt(index, 0);
    var p = 4;
    for (var b = 0; b < nbin; ++b) {
        var bin = readInt(index, p);
        var nchnk = readInt(index, p+4);
//        dlog('bin=' + bin + '; nchnk=' + nchnk);
        p += 8;
        if (intBins[bin]) {
            for (var c = 0; c < nchnk; ++c) {
                var cs = readVob(index, p);
                var ce = readVob(index, p + 8);
                (bin < 4681 ? otherChunks : leafChunks).push(new Chunk(cs, ce));
                p += 16;
            }
        } else {
            p +=  (nchnk * 16);
        }
    }
    // console.log('leafChunks = ' + JSON.stringify(leafChunks));
    // console.log('otherChunks = ' + JSON.stringify(otherChunks));

    var nintv = readInt(index, p);
    // console.log('nintv=' + nintv);
    var lowest = null;
    var minLin = Math.min(min>>14, nintv - 1), maxLin = Math.min(max>>14, nintv - 1);
    for (var i = minLin; i <= maxLin; ++i) {
        var lb =  readVob(index, p + 4 + (i * 8));
        if (!lb) {
            continue;
        }
        if (!lowest || lb.block < lowest.block || (lb.block == lowest.block && lb.offset < lowest.offset)) {
            lowest = lb;
        }
    }
    // console.log('Lowest LB = ' + lowest);

    var prunedOtherChunks = [];
    if (lowest != null) {
        for (var i = 0; i < otherChunks.length; ++i) {
            var chnk = otherChunks[i];
            if (chnk.maxv.block > lowest.block || (chnk.maxv.block == lowest.block && chnk.maxv.offset >= lowest.offset)) {
                prunedOtherChunks.push(chnk);
            }
        }
    }
    // console.log('prunedOtherChunks = ' + JSON.stringify(prunedOtherChunks));
    otherChunks = prunedOtherChunks;

    var intChunks = [];
    for (var i = 0; i < otherChunks.length; ++i) {
        intChunks.push(otherChunks[i]);
    }
    for (var i = 0; i < leafChunks.length; ++i) {
        intChunks.push(leafChunks[i]);
    }

    intChunks.sort(function(c0, c1) {
        var dif = c0.minv.block - c1.minv.block;
        if (dif != 0) {
            return dif;
        } else {
            return c0.minv.offset - c1.minv.offset;
        }
    });
    var mergedChunks = [];
    if (intChunks.length > 0) {
        var cur = intChunks[0];
        for (var i = 1; i < intChunks.length; ++i) {
            var nc = intChunks[i];
            if (nc.minv.block == cur.maxv.block /* && nc.minv.offset == cur.maxv.offset */) { // no point splitting mid-block
                cur = new Chunk(cur.minv, nc.maxv);
            } else {
                mergedChunks.push(cur);
                cur = nc;
            }
        }
        mergedChunks.push(cur);
    }
    // console.log('mergedChunks = ' + JSON.stringify(mergedChunks));

    return mergedChunks;
}

CramFile.prototype.fetch = function(chr, min, max, callback, opts) {
    var thisB = this;
    opts = opts || {};

    var chrId = this.chrToIndex[chr];
    var chunks;
    if (chrId === undefined) {
        chunks = [];
    } else {
        // Fetch this portion of the CRAI if it hasn't been loaded yet.
        if (this.indices[chrId] === null && this.indexChunks.chunks[chrId]) {
            var start_stop = this.indexChunks.chunks[chrId];
            return this.crai.slice(start_stop[0], start_stop[1]).fetch(function(data) {
                var buffer = new Uint8Array(data);
                this.indices[chrId] = buffer;
                return this.fetch(chr, min, max, callback, opts);
            }.bind(this));
        }

        chunks = this.blocksForRange(chrId, min, max);
        if (!chunks) {
            callback(null, 'Error in index fetch');
        }
    }
    
    var records = [];
    var index = 0;
    var data;

    function tramp() {
        if (index >= chunks.length) {
            return callback(records);
        } else if (!data) {
            var c = chunks[index];
            var fetchMin = c.minv.block;
            var fetchMax = c.maxv.block + (1<<16); // *sigh*
            // console.log('fetching ' + fetchMin + ':' + fetchMax);
            thisB.data.slice(fetchMin, fetchMax - fetchMin).fetch(function(r) {
                data = unbgzf(r, c.maxv.block - c.minv.block + 1);
                return tramp();
            });
        } else {
            var ba = new Uint8Array(data);
            var finished = thisB.readCramRecords(ba, chunks[index].minv.offset, records, min, max, chrId, opts);
            data = null;
            ++index;
            if (finished)
                return callback(records);
            else
                return tramp();
        }
    }
    tramp();
}

var SEQRET_DECODER = ['=', 'A', 'C', 'x', 'G', 'x', 'x', 'x', 'T', 'x', 'x', 'x', 'x', 'x', 'x', 'N'];
var CIGAR_DECODER = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', '?', '?', '?', '?', '?', '?', '?'];

function CramRecord() {
}

CramFile.prototype.readCramRecords = function(ba, offset, sink, min, max, chrId, opts) {
    while (true) {
        var blockSize = readInt(ba, offset);
        var blockEnd = offset + blockSize + 4;
        if (blockEnd > ba.length) {
            return false;
        }

        var record = new CramRecord();

        var refID = readInt(ba, offset + 4);
        var pos = readInt(ba, offset + 8);
        
        var bmn = readInt(ba, offset + 12);
        var bin = (bmn & 0xffff0000) >> 16;
        var mq = (bmn & 0xff00) >> 8;
        var nl = bmn & 0xff;

        var flag_nc = readInt(ba, offset + 16);
        var flag = (flag_nc & 0xffff0000) >> 16;
        var nc = flag_nc & 0xffff;
    
        var lseq = readInt(ba, offset + 20);
        
        var nextRef  = readInt(ba, offset + 24);
        var nextPos = readInt(ba, offset + 28);
        
        var tlen = readInt(ba, offset + 32);
    
        record.segment = this.indexToChr[refID];
        record.flag = flag;
        record.pos = pos;
        record.mq = mq;
        if (opts.light)
            record.seqLength = lseq;

        if (!opts.light || opts.includeName) {
            var readName = '';
            for (var j = 0; j < nl-1; ++j) {
                readName += String.fromCharCode(ba[offset + 36 + j]);
            }
            record.readName = readName;
        }
        
        if (!opts.light) {
            if (nextRef >= 0) {
                record.nextSegment = this.indexToChr[nextRef];
                record.nextPos = nextPos;
            }
        
            var p = offset + 36 + nl;

            var cigar = '';
            for (var c = 0; c < nc; ++c) {
                var cigop = readInt(ba, p);
                cigar = cigar + (cigop>>4) + CIGAR_DECODER[cigop & 0xf];
                p += 4;
            }
            record.cigar = cigar;
        
            var seq = '';
            var seqBytes = (lseq + 1) >> 1;
            for (var j = 0; j < seqBytes; ++j) {
                var sb = ba[p + j];
                seq += SEQRET_DECODER[(sb & 0xf0) >> 4];
                if (seq.length < lseq)
                    seq += SEQRET_DECODER[(sb & 0x0f)];
            }
            p += seqBytes;
            record.seq = seq;

            var qseq = '';
            for (var j = 0; j < lseq; ++j) {
                qseq += String.fromCharCode(ba[p + j] + 33);
            }
            p += lseq;
            record.quals = qseq;

            while (p < blockEnd) {
                var tag = String.fromCharCode(ba[p], ba[p + 1]);
                var type = String.fromCharCode(ba[p + 2]);
                var value;

                if (type == 'A') {
                    value = String.fromCharCode(ba[p + 3]);
                    p += 4;
                } else if (type == 'i' || type == 'I') {
                    value = readInt(ba, p + 3);
                    p += 7;
                } else if (type == 'c' || type == 'C') {
                    value = ba[p + 3];
                    p += 4;
                } else if (type == 's' || type == 'S') {
                    value = readShort(ba, p + 3);
                    p += 5;
                } else if (type == 'f') {
                    value = readFloat(ba, p + 3);
                    p += 7;
                } else if (type == 'Z' || type == 'H') {
                    p += 3;
                    value = '';
                    for (;;) {
                        var cc = ba[p++];
                        if (cc == 0) {
                            break;
                        } else {
                            value += String.fromCharCode(cc);
                        }
                    }
                } else if (type == 'B') {
                    var atype = String.fromCharCode(ba[p + 3]);
                    var alen = readInt(ba, p + 4);
                    var elen;
                    var reader;
                    if (atype == 'i' || atype == 'I' || atype == 'f') {
                        elen = 4;
                        if (atype == 'f')
                            reader = readFloat;
                        else
                            reader = readInt;
                    } else if (atype == 's' || atype == 'S') {
                        elen = 2;
                        reader = readShort;
                    } else if (atype == 'c' || atype == 'C') {
                        elen = 1;
                        reader = readByte;
                    } else {
                        throw 'Unknown array type ' + atype;
                    }

                    p += 8;
                    value = [];
                    for (var i = 0; i < alen; ++i) {
                        value.push(reader(ba, p));
                        p += elen;
                    }
                } else {
                    throw 'Unknown type '+ type;
                }
                record[tag] = value;
            }
        }

        if (!min || record.pos <= max && record.pos + lseq >= min) {
            if (chrId === undefined || refID == chrId) {
                sink.push(record);
            }
        }
        if (record.pos > max) {
            return true;
        }
        offset = blockEnd;
    }

    // Exits via top of loop.
};

if (typeof(module) !== 'undefined') {
    module.exports = {
        makeCram: makeCram,
        CRAM_MAGIC: CRAM_MAGIC,
        CRAI_MAGIC: CRAI_MAGIC,
        CramFlags: CramFlags
    };
}
