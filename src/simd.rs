use num_traits::PrimInt;
#[cfg(any(target_arch = "aarch64", target_arch = "x86_64"))]
use std::any::TypeId;

pub(crate) const BLOCK_SIZE: usize = 32;

#[derive(Clone, Copy, Debug)]
pub(crate) enum MaskBackend {
    Scalar,
    #[cfg(target_arch = "aarch64")]
    Neon,
    #[cfg(target_arch = "x86_64")]
    Avx2,
}

#[inline]
pub(crate) fn detect_backend() -> MaskBackend {
    #[cfg(target_arch = "aarch64")]
    {
        return MaskBackend::Neon;
    }

    #[cfg(target_arch = "x86_64")]
    {
        if std::is_x86_feature_detected!("avx2") {
            return MaskBackend::Avx2;
        }
    }

    #[allow(unreachable_code)]
    MaskBackend::Scalar
}

#[inline(always)]
pub(crate) fn overlap_mask<I>(
    backend: MaskBackend,
    starts: &[I],
    stops: &[I],
    query_start: I,
    query_stop: I,
) -> u32
where
    I: PrimInt + 'static,
{
    debug_assert_eq!(starts.len(), stops.len());
    debug_assert!(starts.len() <= BLOCK_SIZE);
    let _ = backend;

    #[cfg(target_arch = "aarch64")]
    if matches!(backend, MaskBackend::Neon) {
        macro_rules! neon_as {
            ($source:ty, $repr:ty, $function:path) => {
                if TypeId::of::<I>() == TypeId::of::<$source>() {
                    // TypeId proves the source type. On AArch64, usize/isize have the
                    // same representation and alignment as u64/i64 respectively.
                    let starts = unsafe {
                        std::slice::from_raw_parts(starts.as_ptr().cast::<$repr>(), starts.len())
                    };
                    let stops = unsafe {
                        std::slice::from_raw_parts(stops.as_ptr().cast::<$repr>(), stops.len())
                    };
                    let query_start = unsafe { *(&query_start as *const I).cast::<$repr>() };
                    let query_stop = unsafe { *(&query_stop as *const I).cast::<$repr>() };
                    return unsafe { $function(starts, stops, query_start, query_stop) };
                }
            };
        }

        neon_as!(u8, u8, neon::mask_u8);
        neon_as!(i8, i8, neon::mask_i8);
        neon_as!(u16, u16, neon::mask_u16);
        neon_as!(i16, i16, neon::mask_i16);
        neon_as!(u32, u32, neon::mask_u32);
        neon_as!(i32, i32, neon::mask_i32);
        neon_as!(u64, u64, neon::mask_u64);
        neon_as!(i64, i64, neon::mask_i64);
        neon_as!(usize, u64, neon::mask_u64);
        neon_as!(isize, i64, neon::mask_i64);
    }

    #[cfg(target_arch = "x86_64")]
    if matches!(backend, MaskBackend::Avx2) {
        macro_rules! avx2_as {
            ($source:ty, $repr:ty, $function:path) => {
                if TypeId::of::<I>() == TypeId::of::<$source>() {
                    // TypeId proves the source type. On x86-64, usize/isize have the
                    // same representation and alignment as u64/i64 respectively.
                    let starts = unsafe {
                        std::slice::from_raw_parts(starts.as_ptr().cast::<$repr>(), starts.len())
                    };
                    let stops = unsafe {
                        std::slice::from_raw_parts(stops.as_ptr().cast::<$repr>(), stops.len())
                    };
                    let query_start = unsafe { *(&query_start as *const I).cast::<$repr>() };
                    let query_stop = unsafe { *(&query_stop as *const I).cast::<$repr>() };
                    return unsafe { $function(starts, stops, query_start, query_stop) };
                }
            };
        }

        avx2_as!(u8, u8, avx2::mask_u8);
        avx2_as!(i8, i8, avx2::mask_i8);
        avx2_as!(u16, u16, avx2::mask_u16);
        avx2_as!(i16, i16, avx2::mask_i16);
        avx2_as!(u32, u32, avx2::mask_u32);
        avx2_as!(i32, i32, avx2::mask_i32);
        avx2_as!(u64, u64, avx2::mask_u64);
        avx2_as!(i64, i64, avx2::mask_i64);
        avx2_as!(usize, u64, avx2::mask_u64);
        avx2_as!(isize, i64, avx2::mask_i64);
    }

    scalar_mask(starts, stops, query_start, query_stop)
}

#[inline(always)]
fn scalar_mask<I: PrimInt>(starts: &[I], stops: &[I], query_start: I, query_stop: I) -> u32 {
    let mut mask = 0_u32;
    for lane in 0..starts.len() {
        if stops[lane] > query_start && starts[lane] < query_stop {
            mask |= 1 << lane;
        }
    }
    mask
}

#[cfg(target_arch = "aarch64")]
mod neon {
    use std::arch::aarch64::*;

    #[inline(always)]
    unsafe fn bits_u8x16(overlap: uint8x16_t) -> u16 {
        let weights = vld1_u8([1_u8, 2, 4, 8, 16, 32, 64, 128].as_ptr());
        let low = vaddv_u8(vand_u8(vget_low_u8(overlap), weights));
        let high = vaddv_u8(vand_u8(vget_high_u8(overlap), weights));
        u16::from(low) | (u16::from(high) << 8)
    }

    #[inline(always)]
    unsafe fn bits_u16x8(overlap: uint16x8_t) -> u8 {
        let weights = vld1_u8([1_u8, 2, 4, 8, 16, 32, 64, 128].as_ptr());
        vaddv_u8(vand_u8(vmovn_u16(overlap), weights))
    }

    #[inline(always)]
    unsafe fn bits_u16x16(low: uint16x8_t, high: uint16x8_t) -> u16 {
        bits_u8x16(vcombine_u8(vmovn_u16(low), vmovn_u16(high)))
    }

    #[inline(always)]
    unsafe fn bits_u32x4(overlap: uint32x4_t) -> u32 {
        let weights = vld1q_u32([1_u32, 2, 4, 8].as_ptr());
        vaddvq_u32(vandq_u32(overlap, weights))
    }

    #[inline(always)]
    unsafe fn bits_u32x8(low: uint32x4_t, high: uint32x4_t) -> u16 {
        let weights = vld1q_u16([1_u16, 2, 4, 8, 16, 32, 64, 128].as_ptr());
        vaddvq_u16(vandq_u16(
            vcombine_u16(vmovn_u32(low), vmovn_u32(high)),
            weights,
        ))
    }

    #[inline(always)]
    unsafe fn bits_u64x4(low: uint64x2_t, high: uint64x2_t) -> u32 {
        let weights = vld1q_u32([1_u32, 2, 4, 8].as_ptr());
        vaddvq_u32(vandq_u32(
            vcombine_u32(vmovn_u64(low), vmovn_u64(high)),
            weights,
        ))
    }

    #[inline(always)]
    unsafe fn scalar_tail<I: Ord + Copy>(
        starts: &[I],
        stops: &[I],
        query_start: I,
        query_stop: I,
        mut lane: usize,
        mut mask: u32,
    ) -> u32 {
        while lane < starts.len() {
            if *stops.get_unchecked(lane) > query_start && *starts.get_unchecked(lane) < query_stop
            {
                mask |= 1 << lane;
            }
            lane += 1;
        }
        mask
    }

    pub(super) unsafe fn mask_u8(
        starts: &[u8],
        stops: &[u8],
        query_start: u8,
        query_stop: u8,
    ) -> u32 {
        let query_start_v = vdupq_n_u8(query_start);
        let query_stop_v = vdupq_n_u8(query_stop);
        let mut lane = 0;
        let mut mask = 0;
        while lane + 32 <= starts.len() {
            let lane_starts = vld1q_u8_x2(starts.as_ptr().add(lane));
            let lane_stops = vld1q_u8_x2(stops.as_ptr().add(lane));
            let low = vandq_u8(
                vcgtq_u8(lane_stops.0, query_start_v),
                vcgtq_u8(query_stop_v, lane_starts.0),
            );
            let high = vandq_u8(
                vcgtq_u8(lane_stops.1, query_start_v),
                vcgtq_u8(query_stop_v, lane_starts.1),
            );
            mask |= u32::from(bits_u8x16(low)) << lane;
            mask |= u32::from(bits_u8x16(high)) << (lane + 16);
            lane += 32;
        }
        while lane + 16 <= starts.len() {
            let lane_starts = vld1q_u8(starts.as_ptr().add(lane));
            let lane_stops = vld1q_u8(stops.as_ptr().add(lane));
            let overlap = vandq_u8(
                vcgtq_u8(lane_stops, query_start_v),
                vcgtq_u8(query_stop_v, lane_starts),
            );
            mask |= u32::from(bits_u8x16(overlap)) << lane;
            lane += 16;
        }
        scalar_tail(starts, stops, query_start, query_stop, lane, mask)
    }

    pub(super) unsafe fn mask_i8(
        starts: &[i8],
        stops: &[i8],
        query_start: i8,
        query_stop: i8,
    ) -> u32 {
        let query_start_v = vdupq_n_s8(query_start);
        let query_stop_v = vdupq_n_s8(query_stop);
        let mut lane = 0;
        let mut mask = 0;
        while lane + 32 <= starts.len() {
            let lane_starts = vld1q_s8_x2(starts.as_ptr().add(lane));
            let lane_stops = vld1q_s8_x2(stops.as_ptr().add(lane));
            let low = vandq_u8(
                vcgtq_s8(lane_stops.0, query_start_v),
                vcgtq_s8(query_stop_v, lane_starts.0),
            );
            let high = vandq_u8(
                vcgtq_s8(lane_stops.1, query_start_v),
                vcgtq_s8(query_stop_v, lane_starts.1),
            );
            mask |= u32::from(bits_u8x16(low)) << lane;
            mask |= u32::from(bits_u8x16(high)) << (lane + 16);
            lane += 32;
        }
        while lane + 16 <= starts.len() {
            let lane_starts = vld1q_s8(starts.as_ptr().add(lane));
            let lane_stops = vld1q_s8(stops.as_ptr().add(lane));
            let overlap = vandq_u8(
                vcgtq_s8(lane_stops, query_start_v),
                vcgtq_s8(query_stop_v, lane_starts),
            );
            mask |= u32::from(bits_u8x16(overlap)) << lane;
            lane += 16;
        }
        scalar_tail(starts, stops, query_start, query_stop, lane, mask)
    }

    pub(super) unsafe fn mask_u16(
        starts: &[u16],
        stops: &[u16],
        query_start: u16,
        query_stop: u16,
    ) -> u32 {
        let query_start_v = vdupq_n_u16(query_start);
        let query_stop_v = vdupq_n_u16(query_stop);
        let mut lane = 0;
        let mut mask = 0;
        while lane + 16 <= starts.len() {
            let lane_starts = vld1q_u16_x2(starts.as_ptr().add(lane));
            let lane_stops = vld1q_u16_x2(stops.as_ptr().add(lane));
            let low = vandq_u16(
                vcgtq_u16(lane_stops.0, query_start_v),
                vcgtq_u16(query_stop_v, lane_starts.0),
            );
            let high = vandq_u16(
                vcgtq_u16(lane_stops.1, query_start_v),
                vcgtq_u16(query_stop_v, lane_starts.1),
            );
            mask |= u32::from(bits_u16x16(low, high)) << lane;
            lane += 16;
        }
        while lane + 8 <= starts.len() {
            let lane_starts = vld1q_u16(starts.as_ptr().add(lane));
            let lane_stops = vld1q_u16(stops.as_ptr().add(lane));
            let overlap = vandq_u16(
                vcgtq_u16(lane_stops, query_start_v),
                vcgtq_u16(query_stop_v, lane_starts),
            );
            mask |= u32::from(bits_u16x8(overlap)) << lane;
            lane += 8;
        }
        scalar_tail(starts, stops, query_start, query_stop, lane, mask)
    }

    pub(super) unsafe fn mask_i16(
        starts: &[i16],
        stops: &[i16],
        query_start: i16,
        query_stop: i16,
    ) -> u32 {
        let query_start_v = vdupq_n_s16(query_start);
        let query_stop_v = vdupq_n_s16(query_stop);
        let mut lane = 0;
        let mut mask = 0;
        while lane + 16 <= starts.len() {
            let lane_starts = vld1q_s16_x2(starts.as_ptr().add(lane));
            let lane_stops = vld1q_s16_x2(stops.as_ptr().add(lane));
            let low = vandq_u16(
                vcgtq_s16(lane_stops.0, query_start_v),
                vcgtq_s16(query_stop_v, lane_starts.0),
            );
            let high = vandq_u16(
                vcgtq_s16(lane_stops.1, query_start_v),
                vcgtq_s16(query_stop_v, lane_starts.1),
            );
            mask |= u32::from(bits_u16x16(low, high)) << lane;
            lane += 16;
        }
        while lane + 8 <= starts.len() {
            let lane_starts = vld1q_s16(starts.as_ptr().add(lane));
            let lane_stops = vld1q_s16(stops.as_ptr().add(lane));
            let overlap = vandq_u16(
                vcgtq_s16(lane_stops, query_start_v),
                vcgtq_s16(query_stop_v, lane_starts),
            );
            mask |= u32::from(bits_u16x8(overlap)) << lane;
            lane += 8;
        }
        scalar_tail(starts, stops, query_start, query_stop, lane, mask)
    }

    pub(super) unsafe fn mask_u32(
        starts: &[u32],
        stops: &[u32],
        query_start: u32,
        query_stop: u32,
    ) -> u32 {
        let query_start_v = vdupq_n_u32(query_start);
        let query_stop_v = vdupq_n_u32(query_stop);
        let simd_len = starts.len() / 4 * 4;
        let mut lane = 0;
        let mut mask = 0;
        while lane + 8 <= simd_len {
            let lane_starts = vld1q_u32_x2(starts.as_ptr().add(lane));
            let lane_stops = vld1q_u32_x2(stops.as_ptr().add(lane));
            let low = vandq_u32(
                vcgtq_u32(lane_stops.0, query_start_v),
                vcgtq_u32(query_stop_v, lane_starts.0),
            );
            let high = vandq_u32(
                vcgtq_u32(lane_stops.1, query_start_v),
                vcgtq_u32(query_stop_v, lane_starts.1),
            );
            mask |= u32::from(bits_u32x8(low, high)) << lane;
            lane += 8;
        }
        while lane < simd_len {
            let lane_starts = vld1q_u32(starts.as_ptr().add(lane));
            let lane_stops = vld1q_u32(stops.as_ptr().add(lane));
            let overlap = vandq_u32(
                vcgtq_u32(lane_stops, query_start_v),
                vcgtq_u32(query_stop_v, lane_starts),
            );
            mask |= bits_u32x4(overlap) << lane;
            lane += 4;
        }
        scalar_tail(starts, stops, query_start, query_stop, lane, mask)
    }

    pub(super) unsafe fn mask_i32(
        starts: &[i32],
        stops: &[i32],
        query_start: i32,
        query_stop: i32,
    ) -> u32 {
        let query_start_v = vdupq_n_s32(query_start);
        let query_stop_v = vdupq_n_s32(query_stop);
        let simd_len = starts.len() / 4 * 4;
        let mut lane = 0;
        let mut mask = 0;
        while lane + 8 <= simd_len {
            let lane_starts = vld1q_s32_x2(starts.as_ptr().add(lane));
            let lane_stops = vld1q_s32_x2(stops.as_ptr().add(lane));
            let low = vandq_u32(
                vcgtq_s32(lane_stops.0, query_start_v),
                vcgtq_s32(query_stop_v, lane_starts.0),
            );
            let high = vandq_u32(
                vcgtq_s32(lane_stops.1, query_start_v),
                vcgtq_s32(query_stop_v, lane_starts.1),
            );
            mask |= u32::from(bits_u32x8(low, high)) << lane;
            lane += 8;
        }
        while lane < simd_len {
            let lane_starts = vld1q_s32(starts.as_ptr().add(lane));
            let lane_stops = vld1q_s32(stops.as_ptr().add(lane));
            let overlap = vandq_u32(
                vcgtq_s32(lane_stops, query_start_v),
                vcgtq_s32(query_stop_v, lane_starts),
            );
            mask |= bits_u32x4(overlap) << lane;
            lane += 4;
        }
        scalar_tail(starts, stops, query_start, query_stop, lane, mask)
    }

    pub(super) unsafe fn mask_u64(
        starts: &[u64],
        stops: &[u64],
        query_start: u64,
        query_stop: u64,
    ) -> u32 {
        let query_start_v = vdupq_n_u64(query_start);
        let query_stop_v = vdupq_n_u64(query_stop);
        let mut lane = 0;
        let mut mask = 0;
        while lane + 4 <= starts.len() {
            let lane_starts = vld1q_u64_x2(starts.as_ptr().add(lane));
            let lane_stops = vld1q_u64_x2(stops.as_ptr().add(lane));
            let low = vandq_u64(
                vcgtq_u64(lane_stops.0, query_start_v),
                vcgtq_u64(query_stop_v, lane_starts.0),
            );
            let high = vandq_u64(
                vcgtq_u64(lane_stops.1, query_start_v),
                vcgtq_u64(query_stop_v, lane_starts.1),
            );
            mask |= bits_u64x4(low, high) << lane;
            lane += 4;
        }
        scalar_tail(starts, stops, query_start, query_stop, lane, mask)
    }

    pub(super) unsafe fn mask_i64(
        starts: &[i64],
        stops: &[i64],
        query_start: i64,
        query_stop: i64,
    ) -> u32 {
        let query_start_v = vdupq_n_s64(query_start);
        let query_stop_v = vdupq_n_s64(query_stop);
        let mut lane = 0;
        let mut mask = 0;
        while lane + 4 <= starts.len() {
            let lane_starts = vld1q_s64_x2(starts.as_ptr().add(lane));
            let lane_stops = vld1q_s64_x2(stops.as_ptr().add(lane));
            let low = vandq_u64(
                vcgtq_s64(lane_stops.0, query_start_v),
                vcgtq_s64(query_stop_v, lane_starts.0),
            );
            let high = vandq_u64(
                vcgtq_s64(lane_stops.1, query_start_v),
                vcgtq_s64(query_stop_v, lane_starts.1),
            );
            mask |= bits_u64x4(low, high) << lane;
            lane += 4;
        }
        scalar_tail(starts, stops, query_start, query_stop, lane, mask)
    }
}

#[cfg(target_arch = "x86_64")]
mod avx2 {
    use std::arch::x86_64::*;

    #[inline(always)]
    fn compact_u16_movemask(mut bits: u32) -> u32 {
        bits &= 0x5555_5555;
        bits = (bits | (bits >> 1)) & 0x3333_3333;
        bits = (bits | (bits >> 2)) & 0x0f0f_0f0f;
        bits = (bits | (bits >> 4)) & 0x00ff_00ff;
        (bits | (bits >> 8)) & 0x0000_ffff
    }

    #[inline(always)]
    unsafe fn scalar_tail<I: Ord + Copy>(
        starts: &[I],
        stops: &[I],
        query_start: I,
        query_stop: I,
        mut lane: usize,
        mut mask: u32,
    ) -> u32 {
        while lane < starts.len() {
            if *stops.get_unchecked(lane) > query_start && *starts.get_unchecked(lane) < query_stop
            {
                mask |= 1 << lane;
            }
            lane += 1;
        }
        mask
    }

    #[target_feature(enable = "avx2")]
    pub(super) unsafe fn mask_u8(
        starts: &[u8],
        stops: &[u8],
        query_start: u8,
        query_stop: u8,
    ) -> u32 {
        let bias = _mm256_set1_epi8(i8::MIN);
        let query_start_v = _mm256_xor_si256(_mm256_set1_epi8(query_start as i8), bias);
        let query_stop_v = _mm256_xor_si256(_mm256_set1_epi8(query_stop as i8), bias);
        let mut lane = 0;
        let mut mask = 0;
        while lane + 32 <= starts.len() {
            let lane_starts = _mm256_loadu_si256(starts.as_ptr().add(lane).cast());
            let lane_stops = _mm256_loadu_si256(stops.as_ptr().add(lane).cast());
            let overlap = _mm256_and_si256(
                _mm256_cmpgt_epi8(_mm256_xor_si256(lane_stops, bias), query_start_v),
                _mm256_cmpgt_epi8(query_stop_v, _mm256_xor_si256(lane_starts, bias)),
            );
            mask |= (_mm256_movemask_epi8(overlap) as u32) << lane;
            lane += 32;
        }
        scalar_tail(starts, stops, query_start, query_stop, lane, mask)
    }

    #[target_feature(enable = "avx2")]
    pub(super) unsafe fn mask_i8(
        starts: &[i8],
        stops: &[i8],
        query_start: i8,
        query_stop: i8,
    ) -> u32 {
        let query_start_v = _mm256_set1_epi8(query_start);
        let query_stop_v = _mm256_set1_epi8(query_stop);
        let mut lane = 0;
        let mut mask = 0;
        while lane + 32 <= starts.len() {
            let lane_starts = _mm256_loadu_si256(starts.as_ptr().add(lane).cast());
            let lane_stops = _mm256_loadu_si256(stops.as_ptr().add(lane).cast());
            let overlap = _mm256_and_si256(
                _mm256_cmpgt_epi8(lane_stops, query_start_v),
                _mm256_cmpgt_epi8(query_stop_v, lane_starts),
            );
            mask |= (_mm256_movemask_epi8(overlap) as u32) << lane;
            lane += 32;
        }
        scalar_tail(starts, stops, query_start, query_stop, lane, mask)
    }

    #[target_feature(enable = "avx2")]
    pub(super) unsafe fn mask_u16(
        starts: &[u16],
        stops: &[u16],
        query_start: u16,
        query_stop: u16,
    ) -> u32 {
        let bias = _mm256_set1_epi16(i16::MIN);
        let query_start_v = _mm256_xor_si256(_mm256_set1_epi16(query_start as i16), bias);
        let query_stop_v = _mm256_xor_si256(_mm256_set1_epi16(query_stop as i16), bias);
        let mut lane = 0;
        let mut mask = 0;
        while lane + 16 <= starts.len() {
            let lane_starts = _mm256_loadu_si256(starts.as_ptr().add(lane).cast());
            let lane_stops = _mm256_loadu_si256(stops.as_ptr().add(lane).cast());
            let overlap = _mm256_and_si256(
                _mm256_cmpgt_epi16(_mm256_xor_si256(lane_stops, bias), query_start_v),
                _mm256_cmpgt_epi16(query_stop_v, _mm256_xor_si256(lane_starts, bias)),
            );
            mask |= compact_u16_movemask(_mm256_movemask_epi8(overlap) as u32) << lane;
            lane += 16;
        }
        scalar_tail(starts, stops, query_start, query_stop, lane, mask)
    }

    #[target_feature(enable = "avx2")]
    pub(super) unsafe fn mask_i16(
        starts: &[i16],
        stops: &[i16],
        query_start: i16,
        query_stop: i16,
    ) -> u32 {
        let query_start_v = _mm256_set1_epi16(query_start);
        let query_stop_v = _mm256_set1_epi16(query_stop);
        let mut lane = 0;
        let mut mask = 0;
        while lane + 16 <= starts.len() {
            let lane_starts = _mm256_loadu_si256(starts.as_ptr().add(lane).cast());
            let lane_stops = _mm256_loadu_si256(stops.as_ptr().add(lane).cast());
            let overlap = _mm256_and_si256(
                _mm256_cmpgt_epi16(lane_stops, query_start_v),
                _mm256_cmpgt_epi16(query_stop_v, lane_starts),
            );
            mask |= compact_u16_movemask(_mm256_movemask_epi8(overlap) as u32) << lane;
            lane += 16;
        }
        scalar_tail(starts, stops, query_start, query_stop, lane, mask)
    }

    #[target_feature(enable = "avx2")]
    pub(super) unsafe fn mask_u32(
        starts: &[u32],
        stops: &[u32],
        query_start: u32,
        query_stop: u32,
    ) -> u32 {
        let bias = _mm256_set1_epi32(i32::MIN);
        let query_start_v = _mm256_xor_si256(_mm256_set1_epi32(query_start as i32), bias);
        let query_stop_v = _mm256_xor_si256(_mm256_set1_epi32(query_stop as i32), bias);
        let mut lane = 0;
        let mut mask = 0;
        while lane + 8 <= starts.len() {
            let lane_starts = _mm256_loadu_si256(starts.as_ptr().add(lane).cast());
            let lane_stops = _mm256_loadu_si256(stops.as_ptr().add(lane).cast());
            let overlap = _mm256_and_si256(
                _mm256_cmpgt_epi32(_mm256_xor_si256(lane_stops, bias), query_start_v),
                _mm256_cmpgt_epi32(query_stop_v, _mm256_xor_si256(lane_starts, bias)),
            );
            mask |= (_mm256_movemask_ps(_mm256_castsi256_ps(overlap)) as u32) << lane;
            lane += 8;
        }
        scalar_tail(starts, stops, query_start, query_stop, lane, mask)
    }

    #[target_feature(enable = "avx2")]
    pub(super) unsafe fn mask_i32(
        starts: &[i32],
        stops: &[i32],
        query_start: i32,
        query_stop: i32,
    ) -> u32 {
        let query_start_v = _mm256_set1_epi32(query_start);
        let query_stop_v = _mm256_set1_epi32(query_stop);
        let mut lane = 0;
        let mut mask = 0;
        while lane + 8 <= starts.len() {
            let lane_starts = _mm256_loadu_si256(starts.as_ptr().add(lane).cast());
            let lane_stops = _mm256_loadu_si256(stops.as_ptr().add(lane).cast());
            let overlap = _mm256_and_si256(
                _mm256_cmpgt_epi32(lane_stops, query_start_v),
                _mm256_cmpgt_epi32(query_stop_v, lane_starts),
            );
            mask |= (_mm256_movemask_ps(_mm256_castsi256_ps(overlap)) as u32) << lane;
            lane += 8;
        }
        scalar_tail(starts, stops, query_start, query_stop, lane, mask)
    }

    #[target_feature(enable = "avx2")]
    pub(super) unsafe fn mask_u64(
        starts: &[u64],
        stops: &[u64],
        query_start: u64,
        query_stop: u64,
    ) -> u32 {
        let bias = _mm256_set1_epi64x(i64::MIN);
        let query_start_v = _mm256_xor_si256(_mm256_set1_epi64x(query_start as i64), bias);
        let query_stop_v = _mm256_xor_si256(_mm256_set1_epi64x(query_stop as i64), bias);
        let mut lane = 0;
        let mut mask = 0;
        while lane + 4 <= starts.len() {
            let lane_starts = _mm256_loadu_si256(starts.as_ptr().add(lane).cast());
            let lane_stops = _mm256_loadu_si256(stops.as_ptr().add(lane).cast());
            let overlap = _mm256_and_si256(
                _mm256_cmpgt_epi64(_mm256_xor_si256(lane_stops, bias), query_start_v),
                _mm256_cmpgt_epi64(query_stop_v, _mm256_xor_si256(lane_starts, bias)),
            );
            mask |= (_mm256_movemask_pd(_mm256_castsi256_pd(overlap)) as u32) << lane;
            lane += 4;
        }
        scalar_tail(starts, stops, query_start, query_stop, lane, mask)
    }

    #[target_feature(enable = "avx2")]
    pub(super) unsafe fn mask_i64(
        starts: &[i64],
        stops: &[i64],
        query_start: i64,
        query_stop: i64,
    ) -> u32 {
        let query_start_v = _mm256_set1_epi64x(query_start);
        let query_stop_v = _mm256_set1_epi64x(query_stop);
        let mut lane = 0;
        let mut mask = 0;
        while lane + 4 <= starts.len() {
            let lane_starts = _mm256_loadu_si256(starts.as_ptr().add(lane).cast());
            let lane_stops = _mm256_loadu_si256(stops.as_ptr().add(lane).cast());
            let overlap = _mm256_and_si256(
                _mm256_cmpgt_epi64(lane_stops, query_start_v),
                _mm256_cmpgt_epi64(query_stop_v, lane_starts),
            );
            mask |= (_mm256_movemask_pd(_mm256_castsi256_pd(overlap)) as u32) << lane;
            lane += 4;
        }
        scalar_tail(starts, stops, query_start, query_stop, lane, mask)
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn u16_movemask_compaction_preserves_every_lane_bit() {
            for expected in 0_u32..=u16::MAX.into() {
                let mut expanded = 0_u32;
                for lane in 0..16 {
                    if expected & (1 << lane) != 0 {
                        expanded |= 0b11 << (lane * 2);
                    }
                }
                assert_eq!(compact_u16_movemask(expanded), expected);
            }
        }

        #[test]
        fn primitive_masks_match_scalar_when_avx2_is_available() {
            if !std::is_x86_feature_detected!("avx2") {
                return;
            }

            macro_rules! check_unsigned {
                ($type:ty, $function:path) => {{
                    let starts: Vec<$type> = (0..32).map(|lane| (lane * 4) as $type).collect();
                    let stops: Vec<$type> = starts.iter().map(|start| *start + 20).collect();
                    let expected =
                        super::super::scalar_mask(&starts, &stops, 37 as $type, 91 as $type);
                    assert_eq!(
                        unsafe { $function(&starts, &stops, 37 as $type, 91 as $type) },
                        expected
                    );
                }};
            }

            macro_rules! check_signed {
                ($type:ty, $function:path) => {{
                    let starts: Vec<$type> = (0..32).map(|lane| (lane * 4 - 64) as $type).collect();
                    let stops: Vec<$type> = starts.iter().map(|start| *start + 20).collect();
                    let expected =
                        super::super::scalar_mask(&starts, &stops, -11 as $type, 43 as $type);
                    assert_eq!(
                        unsafe { $function(&starts, &stops, -11 as $type, 43 as $type) },
                        expected
                    );
                }};
            }

            check_unsigned!(u8, mask_u8);
            check_unsigned!(u16, mask_u16);
            check_unsigned!(u32, mask_u32);
            check_unsigned!(u64, mask_u64);
            check_signed!(i8, mask_i8);
            check_signed!(i16, mask_i16);
            check_signed!(i32, mask_i32);
            check_signed!(i64, mask_i64);
        }
    }
}
