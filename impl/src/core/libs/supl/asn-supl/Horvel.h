/*
 * SPDX-FileCopyrightText: (c) 2003, 2004 Lev Walkin <vlm@lionet.info>. All rights reserved.
 * SPDX-License-Identifier: BSD-1-Clause
 * Generated by asn1c-0.9.22 (http://lionet.info/asn1c)
 * From ASN.1 module "ULP-Components"
 *     found in "../supl-common.asn"
 */

#ifndef _Horvel_H
#define _Horvel_H

#include <asn_application.h>

/* Including external dependencies */
#include <BIT_STRING.h>
#include <constr_SEQUENCE.h>

#ifdef __cplusplus
extern "C"
{
#endif

    /* Horvel */
    typedef struct Horvel
    {
        BIT_STRING_t bearing;
        BIT_STRING_t horspeed;
        /*
         * This type is extensible,
         * possible extensions are below.
         */

        /* Context for parsing across buffer boundaries */
        asn_struct_ctx_t _asn_ctx;
    } Horvel_t;

    /* Implementation */
    extern asn_TYPE_descriptor_t asn_DEF_Horvel;

#ifdef __cplusplus
}
#endif

#endif /* _Horvel_H_ */
#include <asn_internal.h>