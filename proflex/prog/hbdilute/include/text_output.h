/*******************************************************************************
 *  MSU ProFlex, formerly called FIRST, is a software developed to predict and *
 *  analyze protein flexibility. * This source file is a part of MSU ProFlex. *
 *                                                                              *
 *  Copyright (C) 1997 - 2006, Michigan State University. *
 *                                                                              *
 *  This program is free software; you can redistribute to academic users only,
 ** it and/or modify it under the terms of the GNU General Public License, *
 *  version 2, as published by the Free Software Foundation. *
 *                                                                              *
 *  This program is distributed in the hope that it will be useful, * but
 *WITHOUT ANY WARRANTY; without even the implied warranty of              *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the * GNU General
 *Public License for more details.                                *
 *                                                                              *
 *  You should have received a copy of the GNU General Public License * along
 *with this program; if not, write to the Free Software                 *
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA, *
 *  or see http://www.gnu.org/licenses/gpl.txt *
 *******************************************************************************/

/*
 * File:   text_output.h
 * Author: Sameer Arora
 *
 * Created on November 10, 2003, 5:23 PM
 */

#ifndef _text_output_H
#define _text_output_H

#include "types.h"
#include <stdio.h>

void print_text_decomp(clusters *[], int, int, FILE *, int[], int, int, float,
                       float);

#endif /* _text_output_H */
