import { NextResponse } from 'next/server'
import { query } from '@/lib/db'

export async function GET() {
  try {
    // Get total visit count
    const visitCountResult = await query(
      'SELECT COALESCE(SUM(visit_count), 0) as total FROM visits'
    )
    const visitCount = parseInt(visitCountResult.rows[0].total)

    // Get unique visitors count
    const uniqueVisitorsResult = await query('SELECT COUNT(*) as count FROM visits')
    const uniqueVisitors = parseInt(uniqueVisitorsResult.rows[0].count)

    return NextResponse.json({
      visitCount,
      uniqueVisitors,
    })
  } catch (error) {
    console.error('Error fetching visit stats:', error)
    return NextResponse.json(
      { visitCount: 0, uniqueVisitors: 0 },
      { status: 500 }
    )
  }
}
